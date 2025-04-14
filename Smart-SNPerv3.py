F = input('Please select the SNP design method{A (enter the wild-type sequence and mutation sequence) B (enter the RS number and crawl the network database for design) C (enter the RS number and utilize the local database for design)}:')
if F == 'A':
    import pandas as pd
    import re
    import os
    import RNA


    def reverse_complement(sequence):
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return ''.join(complement[base] for base in reversed(sequence.upper()))


    def read_fasta_sequence(file_names):
        sequence = []
        for file_name in file_names:
            with open(file_name, 'r') as file:
                for line in file:
                    if not line.startswith('>'):
                        cleaned_line = re.sub(r'[^ATCGN]', '', line.strip())
                        sequence.append(cleaned_line)
        return ''.join(sequence)


    def find_differences(str1, str2):
        min_len = min(len(str1), len(str2))
        differences = []
        for i in range(min_len):
            if str1[i].upper() != str2[i].upper():
                differences.append(i)
        if len(str1) != len(str2):
            differences.extend(range(min_len, max(len(str1), len(str2))))
        return differences


    def mutate_string(g, pos):
        mutation_dict = {
            'A': ['C', 'T', 'G'],
            'T': ['C', 'A', 'G'],
            'G': ['C', 'A', 'T'],
            'C': ['A', 'G', 'T']
        }
        g_upper = g.upper()
        pos -= 1
        if len(g_upper) < 5:
            return [g_upper]
        result = []
        for i in range(2, 5):
            if i == pos:
                continue
            char = g_upper[i]
            mutations = mutation_dict.get(char, [char])
            for mutation in mutations:
                mutated_g = g_upper[:i] + mutation + g_upper[i + 1:]
                result.append(mutated_g)
        return result


    def find_pattern_and_report(seq1, seq2, pos, direction):
        results = []
        seq1 = seq1.upper()
        seq2 = seq2.upper()


        stand_symbol = '+' if direction == '正向结果' else '-'

        if seq1[pos] != 'T' and seq2[pos] == 'T':
            if pos < len(seq1) - 3 and seq1[pos + 1:pos + 3] in ["TA", "TC", "TG", "AT", "CT", "GT"]:
                pam = seq2[pos:pos + 4]
                crRNA_sequence = seq2[pos + 4:pos + 22]
                results.append({
                    'stand': stand_symbol,
                    'pos in sequence': pos + 1,
                    'PAM': pam,
                    'protospacer': crRNA_sequence,
                    'pos in PAM': 1,
                })
            if pos > 1 and seq1[pos - 2:pos] in ["TA", "TC", "TG", "AT", "CT", "GT"]:
                pam = seq2[pos - 2:pos + 2]
                crRNA_sequence = seq2[pos + 2:pos + 20]
                results.append({
                    'stand': stand_symbol,
                    'pos in sequence': pos + 1,
                    'PAM': pam,
                    'protospacer': crRNA_sequence,
                    'pos in PAM': 3,
                })
            if pos > 0 and pos < len(seq1) - 2 and seq1[pos - 1] == 'T' and seq1[pos + 1] in ["A", "C", "G"]:
                pam = seq2[pos - 1:pos + 3]
                crRNA_sequence = seq2[pos + 3:pos + 21]
                results.append({
                    'stand': stand_symbol,
                    'pos in sequence': pos + 1,
                    'PAM': pam,
                    'protospacer': crRNA_sequence,
                    'pos in PAM': 2,
                })
            if pos < len(seq1) - 1 and seq1[pos + 1] == 'T' and seq1[pos - 1] in ["A", "C", "G"]:
                pam = seq2[pos - 1:pos + 3]
                crRNA_sequence = seq2[pos + 3:pos + 21]
                results.append({
                    'stand': stand_symbol,
                    'pos in sequence': pos + 1,
                    'PAM': pam,
                    'protospacer': crRNA_sequence,
                    'pos in PAM': 2,
                })

        patterns = [r'ttt[atcg]', r't[atcg]t[atcg]', r'[atcg]tt[atcg]', r'tt[atcg][atcg]']
        start = max(0, pos - 14)
        context = seq2[start:pos]

        for pattern in patterns:
            i = 0
            while i <= len(context) - 4:
                match = re.match(pattern, context[i:], re.IGNORECASE)
                if match:
                    end_pos = match.end() + i + start
                    following_chars = seq2[end_pos:end_pos + 18]
                    new_pos = pos - end_pos + 1
                    result = mutate_string(following_chars, new_pos)
                    results.append({
                        'stand': stand_symbol,
                        'pos in sequence': pos + 1,
                        'PAM': match.group().upper(),
                        'protospacer': following_chars,
                        'pos in protospacer': new_pos,
                        'guide sequence': result
                    })
                    i += 1
                else:
                    i += 1

        return results


    def process_sequences(seq1, seq2, direction):
        differences = find_differences(seq1, seq2)
        results = []
        for pos in differences:
            results.extend(find_pattern_and_report(seq1, seq2, pos, direction))
        return results


    def find_pam_sites(seq1):
        pam_patterns = ["TTTN", "NTTN", "TNTN", "TTNN"]
        pam_sites = []

        seq1 = seq1.upper()

        for i in range(len(seq1) - 3):
            pam_candidate = seq1[i:i + 4]
            for pattern in pam_patterns:
                if all(p == 'N' or p == c for p, c in zip(pattern, pam_candidate)):
                    pam_sites.append(i)
                    break

        return pam_sites


    def compare_crrna_with_targets(crrna, targets):
        results = []

        crrna = crrna.upper()

        for target in targets:
            target = target.upper()
            match_count = sum(1 for a, b in zip(crrna, target) if a == b)
            match_count_first_five = sum(1 for a, b in zip(crrna[:5], target[:5]) if a == b)
            results.append((match_count, match_count_first_five))
        return results


    def perform_specificity_analysis(crrna, seq1):
        pam_sites = find_pam_sites(seq1)
        target_length = len(crrna) + 4


        targets = [
            (seq1[i:i + 4], seq1[i + 4:i + 22])
            for i in pam_sites
            if i + 22 <= len(seq1) and len(seq1[i:i + 22]) == target_length
        ]


        target_sequences = [pam + target for pam, target in targets]


        results = compare_crrna_with_targets(crrna, target_sequences)

        similarity_counts = {}
        for match_count, match_count_first_five in results:
            key = (match_count_first_five, match_count)
            if key not in similarity_counts:
                similarity_counts[key] = 0
            similarity_counts[key] += 1

        detailed_results = []
        for i, ((pam, target), (match_count, match_count_first_five)) in enumerate(zip(targets, results)):
            detailed_results.append({
                'Index': i + 1,
                'PAM': pam,
                'Target Sequence': target,
                'Total Match Count': match_count,
                'First Five Match Count': match_count_first_five
            })

        summary_results = []
        for (first_five, total), count in similarity_counts.items():
            summary_results.append({
                'First Five Match Count': first_five,
                'Total Match Count': total,
                'Count': count
            })


        summary_results = sorted(summary_results, key=lambda x: (-x['First Five Match Count'], -x['Total Match Count']))

        return detailed_results, summary_results


    def predict_rna_structure(sequence):
        fc = RNA.fold_compound(sequence)
        structure, mfe = fc.mfe()
        return structure, mfe


    def predict_rna_structure_for_guide_sequences(guide_sequences):
        rna_structure_data = []
        prefix_sequence = "TAATTTCTACTAAGTGTAGAT".replace("T", "U")


        if isinstance(guide_sequences, str):
            guide_sequences = [guide_sequences]


        for rna_sequence in guide_sequences:
            rna_sequence = rna_sequence.replace("T", "U")
            full_sequence = prefix_sequence + rna_sequence
            structure, mfe = predict_rna_structure(full_sequence)

            rna_structure_data.append({
                'Full RNA Sequence': full_sequence,
                'Predicted Structure': structure,
                'MFE (kcal/mol)': mfe
            })

        return rna_structure_data


    def get_hamming(seq1, seq2):
        if len(seq1) != len(seq2):
            raise ValueError("Sequences must be of equal length")
        return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


    def generate_primers(seq1, site, length=30):
        up_primers = []
        down_primers = []

        for i in range(max(0, site - 125), max(0, site - 50) - length + 1):
            up_primers.append((seq1[i:i + length].upper(), i))

        for i in range(site + 50, site + 125 - length):
            down_seq = seq1[i:i + length]
            down_primers.append((reverse_complement(down_seq).upper(), i + length - 1))

        return up_primers, down_primers


    def calculate_gc_content(seq):
        seq = seq.upper()
        if len(seq) == 0:
            return 0
        gc_count = seq.count('G') + seq.count('C')
        return gc_count / len(seq)


    def calculate_rna_pairing(primer, rna_seq):
        max_pairing = 0
        for i in range(len(rna_seq) - len(primer) + 1):
            pairing = len(primer) - get_hamming(primer, rna_seq[i:i + len(primer)])
            if pairing > max_pairing:
                max_pairing = pairing
        return max_pairing


    def filter_primers(primers, seq1, rna_seq, threshold=15, max_occurrences=30, gc_min=0.3, gc_max=0.7,
                       max_self_pairing=6, max_rna_pairing=6):
        filtered_primers = []

        for primer, start_pos in primers:
            gc_content = calculate_gc_content(primer)
            if gc_content < gc_min or gc_content > gc_max:
                continue

            reverse_comp = reverse_complement(primer)
            self_pairing_count = len(primer) - get_hamming(primer, reverse_comp)
            if self_pairing_count > max_self_pairing:
                continue

            rna_pairing_count = calculate_rna_pairing(primer, rna_seq)
            if rna_pairing_count > max_rna_pairing:
                continue

            count = 0
            for i in range(len(seq1) - len(primer) + 1):
                if get_hamming(primer, seq1[i:i + len(primer)]) <= (len(primer) - threshold):
                    count += 1
                if count > max_occurrences:
                    break
            if count <= max_occurrences:
                filtered_primers.append((primer, start_pos))

        return filtered_primers


    def filter_primer_pairs(up_primers, down_primers, max_pairing=12):
        filtered_pairs = []

        for up_primer, up_start in up_primers:
            for down_primer, down_end in down_primers:
                pairing_count = len(up_primer) - get_hamming(up_primer, down_primer)
                if pairing_count <= max_pairing:
                    distance = down_end - up_start + 1
                    filtered_pairs.append((up_primer, down_primer, distance))

        return filtered_pairs


    def save_to_excel(filtered_pairs, output_file):
        primer_data = [{'Up-Primer': up, 'Down-Primer': down, 'Distance': distance} for up, down, distance in
                       filtered_pairs]
        df = pd.DataFrame(primer_data)

        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            df.to_excel(writer, sheet_name='Primer Pairs', index=False)

        print(f"Filtered primer pairs saved to {output_file}")


    def save_results_to_excel(results, seq1, seq2, seq5, seq6, base_filename='crRNA_result'):
        if not results:
            print("No crRNA results found.")
            return

        prefix_sequence = "TAATTTCTACTAAGTGTAGAT".replace("T", "U")


        processed_results = set()

        for i, result in enumerate(results):

            unique_id = (result.get('PAM'), result.get('pos in sequence'), result.get('protospacer'))
            if unique_id in processed_results:
                continue

            processed_results.add(unique_id)

            df = pd.DataFrame([result])
            pam = result.get('PAM', 'unknown')
            position = result.get('pos in sequence', i + 1)
            output_file = f'{base_filename}_PAM_{pam}_Pos_{position}_{i + 1}.xlsx'

            try:
                detailed_results, summary_results = perform_specificity_analysis(result['protospacer'], seq5)
                detailed_results_rev, summary_results_rev = perform_specificity_analysis(result['protospacer'], seq6)


                if not detailed_results and not detailed_results_rev:
                    print(f"No valid results for {output_file}. Skipping...")
                    continue


                rna_structure_data = predict_rna_structure_for_guide_sequences(
                    result.get('guide sequence') or result.get('protospacer'))


                site = result['pos in sequence'] - 1
                up_primers, down_primers = generate_primers(seq1, site)
                filtered_up_primers = filter_primers(up_primers, seq5, result['protospacer'])
                filtered_down_primers = filter_primers(down_primers, seq6, result['protospacer'])
                filtered_pairs = filter_primer_pairs(filtered_up_primers, filtered_down_primers)
                primer_df = pd.DataFrame(filtered_pairs, columns=['Up-Primer', 'Down-Primer', 'Distance'])

                with pd.ExcelWriter(output_file) as writer:
                    df.to_excel(writer, sheet_name='crRNA Results', index=False)
                    if detailed_results:
                        pd.DataFrame(detailed_results).to_excel(writer, sheet_name='Detailed Results Forward',
                                                                index=False)
                    if summary_results:
                        pd.DataFrame(summary_results).to_excel(writer, sheet_name='Summary Forward', index=False)
                    if detailed_results_rev:
                        pd.DataFrame(detailed_results_rev).to_excel(writer, sheet_name='Detailed Results Reverse',
                                                                    index=False)
                    if summary_results_rev:
                        pd.DataFrame(summary_results_rev).to_excel(writer, sheet_name='Summary Reverse', index=False)

                    # Save RNA structure prediction
                    pd.DataFrame(rna_structure_data).to_excel(writer, sheet_name='RNA Structure', index=False)

                    # Save primers
                    primer_df.to_excel(writer, sheet_name='Primers', index=False)

                print(f"Result saved to {os.path.abspath(output_file)}")
            except Exception as e:
                print(f"Failed to save result {i + 1}: {e}")


    try:
        sequence1 = input('please input wild type sequence: ')
        sequence2 = input('please input mutant type sequence: ')
        file_names_input = input(
            "Please enter multiple FASTA file names (separated by spaces), or press Enter to use default genomes : ").strip()


        if not file_names_input:
            file_names_input = "chrY.fasta chrX.fasta chrMT.fasta chr1.fasta chr2.fasta chr3.fasta chr4.fasta chr5.fasta chr6.fasta chr7.fasta chr8.fasta chr9.fasta chr10.fasta chr11.fasta chr12.fasta chr13.fasta chr14.fasta chr15.fasta chr16.fasta chr17.fasta chr18.fasta chr19.fasta chr20.fasta chr21.fasta chr22.fasta"


        file_names_list = file_names_input.split()
        sequence5 = read_fasta_sequence(file_names_list)
        sequence6 = reverse_complement(sequence5)
        results_forward = process_sequences(sequence1, sequence2, "正向结果")
        sequence3 = reverse_complement(sequence1)
        sequence4 = reverse_complement(sequence2)
        results_reverse = process_sequences(sequence3, sequence4, "反向结果")

        all_results = results_forward + results_reverse

        save_results_to_excel(all_results, sequence2, sequence4, sequence5, sequence6)

        input("Press Enter to exit...")

    except Exception as e:
        print(f"An error occurred: {e}")
        input("Press Enter to exit...")

if F == 'B':
    from selenium import webdriver
    from selenium.webdriver.common.by import By
    from selenium.webdriver.support.ui import WebDriverWait
    from selenium.webdriver.support import expected_conditions as EC
    from selenium.webdriver.chrome.service import Service
    from webdriver_manager.chrome import ChromeDriverManager
    import time
    import re
    import requests
    from bs4 import BeautifulSoup
    import pandas as pd
    import os
    import RNA


    def fetch_snp_data(rs_id):
        url = f"https://www.ncbi.nlm.nih.gov/snp/{rs_id}"
        headers = {
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/85.0.4183.121 Safari/537.36"
        }

        try:
            response = requests.get(url, headers=headers)
            response.raise_for_status()

            soup = BeautifulSoup(response.content, 'html.parser')


            alleles1 = []
            alleles2 = []
            alleles_tag = soup.find('dt', string='Alleles')
            if alleles_tag:
                alleles_raw = alleles_tag.find_next('dd').get_text(strip=True).replace(' ', '')
                for allele in alleles_raw.split('/'):
                    if '>' in allele:
                        parts = allele.split('>')
                        alleles2.append(parts[0])
                        alleles1.append(parts[1])

            # 合并和验证 alleles2
            if len(set(alleles2)) != 1:
                raise ValueError("alleles2 data error")
            alleles2_merged = alleles2[0]

            return alleles1, alleles2_merged

        except requests.exceptions.HTTPError as http_err:
            print(f"HTTP error occurred: {http_err}")
        except requests.exceptions.RequestException as req_err:
            print(f"Request error occurred: {req_err}")
        except ValueError as val_err:
            print(f"Value error occurred: {val_err}")
        except Exception as e:
            print(f"An unexpected error occurred: {e}")

        return None, None


    def introduce_mutations(seq, alleles, position=201):
        seqs_with_mutations = {}
        for i, allele in enumerate(alleles, 1):
            mutated_seq = list(seq)
            mutated_seq[position - 1] = allele
            seqs_with_mutations[f'seq2.{i}'] = ''.join(mutated_seq)

        return seqs_with_mutations


    def reverse_complement(sequence):
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return ''.join(complement[base] for base in reversed(sequence.upper()))


    def find_differences(str1, str2):
        min_len = min(len(str1), len(str2))
        differences = []
        for i in range(min_len):
            if str1[i].lower() != str2[i].lower():
                differences.append(i)
        if len(str1) != len(str2):
            differences.extend(range(min_len, max(len(str1), len(str2))))
        return differences


    def read_fasta_sequence(file_names):
        sequence = []
        for file_name in file_names:
            with open(file_name, 'r') as file:
                for line in file:
                    if not line.startswith('>'):

                        cleaned_line = re.sub(r'[^ATCG]', '', line.strip())
                        sequence.append(cleaned_line)
        return ''.join(sequence)


    def get_mutation_description(seq1, seq2, differences):
        mutations = []
        for pos in differences:
            mutations.append(f"{seq1[pos]}->{seq2[pos]} at {position_text}")
        return ', '.join(mutations)


    def mutate_string(g, pos):
        mutation_dict = {
            'A': ['C', 'T', 'G'],
            'T': ['C', 'A', 'G'],
            'G': ['C', 'A', 'T'],
            'C': ['A', 'G', 'T']
        }
        original_case = [(ch.isupper(), i) for i, ch in enumerate(g)]
        g_upper = g.upper()
        pos -= 1
        if len(g_upper) < 5:
            return [g]
        result = []
        for i in range(2, 5):
            if i == pos:
                continue
            char = g_upper[i]
            if char in mutation_dict:
                mutations = mutation_dict[char]
            else:
                return [g]
            for mutation in mutations:
                mutated_g = g_upper[:i] + mutation + g_upper[i + 1:]
                mutated_g_with_case = ''.join(
                    (ch.lower() if not ig_upper else ch.upper())
                    for ch, (ig_upper, idx) in zip(mutated_g, original_case)
                )
                result.append(mutated_g_with_case)
        return result


    def find_pattern_and_report(seq1, seq2, pos, direction, mutation_desc):
        results = []
        seq1 = seq1.upper()
        seq2 = seq2.upper()
        stand_symbol = '+' if direction == '正向结果' else '-'
        if seq1[pos] != 'T' and seq2[pos] == 'T':
            if pos < len(seq1) - 3 and seq1[pos + 1:pos + 3] in ["TA", "TC", "TG", "AT", "CT", "GT"]:
                pam = seq2[pos:pos + 4]
                crRNA_sequence = seq2[pos + 4:pos + 22]
                results.append({
                    'stand': stand_symbol,
                    'pos in sequence': pos + 1,
                    'PAM': pam,
                    'protospacer': crRNA_sequence,
                    'pos in PAM': 1,
                    'Description of the mutation': mutation_desc
                })
            if pos > 1 and seq1[pos - 2:pos] in ["TA", "TC", "TG", "AT", "CT", "GT"]:
                pam = seq2[pos - 2:pos + 2]
                crRNA_sequence = seq2[pos + 2:pos + 20]
                results.append({
                    'stand': stand_symbol,
                    'pos in sequence': pos + 1,
                    'PAM': pam,
                    'protospacer': crRNA_sequence,
                    'pos in PAM': 3,
                    'Description of the mutation': mutation_desc
                })
            if pos > 0 and pos < len(seq1) - 2 and seq1[pos - 1] == 'T' and seq1[pos + 1] in ["A", "C", "G"]:
                pam = seq2[pos - 1:pos + 3]
                crRNA_sequence = seq2[pos + 3:pos + 21]
                results.append({
                    'stand': stand_symbol,
                    'pos in sequence': pos + 1,
                    'PAM': pam,
                    'protospacer': crRNA_sequence,
                    'pos in PAM': 2,
                    'Description of the mutation': mutation_desc
                })
            if pos < len(seq1) - 1 and seq1[pos + 1] == 'T' and seq1[pos - 1] in ["A", "C", "G"]:
                pam = seq2[pos - 1:pos + 3]
                crRNA_sequence = seq2[pos + 3:pos + 21]
                results.append({
                    'stand': stand_symbol,
                    'pos in sequence': pos + 1,
                    'PAM': pam,
                    'protospacer': crRNA_sequence,
                    'pos in PAM': 2,
                    'Description of the mutation': mutation_desc
                })

        patterns = [r'ttt[atcg]', r't[atcg]t[atcg]', r'[atcg]tt[atcg]', r'tt[atcg][atcg]']
        start = max(0, pos - 14)
        context = seq2[start:pos]

        for pattern in patterns:
            i = 0
            while i <= len(context) - 4:
                match = re.match(pattern, context[i:], re.IGNORECASE)
                if match:
                    end_pos = match.end() + i + start
                    following_chars = seq2[end_pos:end_pos + 18]
                    new_pos = pos - end_pos + 1
                    result = mutate_string(following_chars, new_pos)
                    results.append({
                        'stand': stand_symbol,
                        'pos in sequence': pos + 1,
                        'PAM': match.group().upper(),
                        'protospacer': following_chars,
                        'pos in protospacer': new_pos,
                        'guide sequence': result,
                        'Description of the mutation': mutation_desc
                    })
                    i += 1
                else:
                    i += 1

        return results


    def process_sequences(seq1, seq2, direction):
        differences = find_differences(seq1, seq2)
        mutation_desc = get_mutation_description(seq1, seq2, differences)
        results = []
        for pos in differences:
            results.extend(find_pattern_and_report(seq1, seq2, pos, direction, mutation_desc))
        return results, mutation_desc


    def find_pam_sites(seq1):
        pam_patterns = ["TTTN", "NTTN", "TNTN", "TTNN"]
        pam_sites = []

        seq1 = seq1.upper()

        for i in range(len(seq1) - 3):
            pam_candidate = seq1[i:i + 4]
            for pattern in pam_patterns:
                if all(p == 'N' or p == c for p, c in zip(pattern, pam_candidate)):
                    pam_sites.append(i)
                    break

        return pam_sites


    def compare_crrna_with_targets(crrna, targets):
        results = []

        crrna = crrna.upper()

        for target in targets:
            target = target.upper()
            match_count = sum(1 for a, b in zip(crrna, target) if a == b)
            match_count_first_five = sum(1 for a, b in zip(crrna[:5], target[:5]) if a == b)
            results.append((match_count, match_count_first_five))
        return results


    def perform_specificity_analysis(crrna, seq1, chunk_size=100000, overlap=30):
        detailed_results = []
        summary_results = {}

        start_time = time.time()


        for start in range(0, len(seq1), chunk_size):
            end = min(start + chunk_size + overlap, len(seq1))
            chunk = seq1[start:end]


            pam_sites = find_pam_sites(chunk)
            targets = [(chunk[i:i + 4], chunk[i + 4:i + 22]) for i in pam_sites if i + 22 <= len(chunk)]

            for pam, target in targets:

                match_count_first_five = sum(1 for a, b in zip(crrna[:5], target[:5]) if a == b)

                if match_count_first_five < 4:
                    continue


                match_count_rest = sum(1 for a, b in zip(crrna[5:], target[5:]) if a == b)


                match_count = match_count_first_five + match_count_rest

                detailed_results.append({
                    'Index': len(detailed_results) + 1,
                    'PAM': pam,
                    'Target Sequence': target,
                    'Total Match Count': match_count,
                    'First Five Match Count': match_count_first_five
                })

                key = (match_count_first_five, match_count)
                if key not in summary_results:
                    summary_results[key] = 0
                summary_results[key] += 1

        summary_results = sorted(
            [{'First Five Match Count': k[0], 'Total Match Count': k[1], 'Count': v} for k, v in
             summary_results.items()],
            key=lambda x: (-x['First Five Match Count'], -x['Total Match Count'])
        )

        end_time = time.time()
        print(f"RNA specificity analysis took {end_time - start_time:.2f} seconds.")

        return detailed_results, summary_results


    def predict_rna_structure(sequence):
        fc = RNA.fold_compound(sequence)
        structure, mfe = fc.mfe()
        return structure, mfe


    def predict_rna_structure_for_guide_sequences(guide_sequences):
        rna_structure_data = []
        prefix_sequence = "TAATTTCTACTAAGTGTAGAT".replace("T", "U")


        if isinstance(guide_sequences, str):
            guide_sequences = [guide_sequences]


        for rna_sequence in guide_sequences:
            rna_sequence = rna_sequence.replace("T", "U")
            full_sequence = prefix_sequence + rna_sequence
            structure, mfe = predict_rna_structure(full_sequence)

            rna_structure_data.append({
                'Full RNA Sequence': full_sequence,
                'Predicted Structure': structure,
                'MFE (kcal/mol)': mfe
            })

        return rna_structure_data


    def get_hamming(seq1, seq2):
        if len(seq1) != len(seq2):
            raise ValueError("Sequences must be of equal length")
        return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


    def generate_primers(seq1, site, length=30):
        up_primers = []
        down_primers = []

        for i in range(max(0, site - 125), max(0, site - 50) - length + 1):
            up_primers.append((seq1[i:i + length].upper(), i))

        for i in range(site + 50, site + 125 - length):
            down_seq = seq1[i:i + length]
            down_primers.append((reverse_complement(down_seq).upper(), i + length - 1))

        return up_primers, down_primers


    def calculate_gc_content(seq):
        gc_count = seq.upper().count('G') + seq.upper().count('C')
        return gc_count / len(seq)


    def calculate_rna_pairing(primer, rna_seq):
        max_pairing = 0
        for i in range(len(rna_seq) - len(primer) + 1):
            pairing = len(primer) - get_hamming(primer, rna_seq[i:i + len(primer)])
            if pairing > max_pairing:
                max_pairing = pairing
        return max_pairing


    def filter_primers(primers, seq1, rev_comp_seq1, rna_seq, threshold=20, max_occurrences=22, gc_min=0.3, gc_max=0.7,
                       max_self_pairing=6, max_rna_pairing=6, seed_length=15, chunk_size=100000, overlap=30):
        filtered_primers = []

        start_time = time.time()

        def find_seed_matches(primer, sequence):
            seed = primer[:seed_length]
            positions = []
            for i in range(len(sequence) - seed_length + 1):
                if sequence[i:i + seed_length] == seed:
                    positions.append(i)
            return positions

        for primer, start_pos in primers:
            gc_content = calculate_gc_content(primer)
            if gc_content < gc_min or gc_content > gc_max:
                continue

            reverse_comp = reverse_complement(primer)
            self_pairing_count = len(primer) - get_hamming(primer, reverse_comp)
            if self_pairing_count > max_self_pairing:
                continue

            rna_pairing_count = calculate_rna_pairing(primer, rna_seq)
            if rna_pairing_count > max_rna_pairing:
                continue

            total_count = 0

            for start in range(0, len(seq1), chunk_size):
                end = min(start + chunk_size + overlap, len(seq1))
                chunk = seq1[start:end]
                rev_comp_chunk = rev_comp_seq1[start:end]

                count = 0
                seed_matches = find_seed_matches(primer, chunk)
                for match_pos in seed_matches:
                    if match_pos + len(primer) <= len(chunk):
                        full_match = chunk[match_pos:match_pos + len(primer)]
                        if get_hamming(primer, full_match) <= (len(primer) - threshold):
                            count += 1
                            total_count += 1
                            if total_count > max_occurrences:
                                break

                if total_count > max_occurrences:
                    break

                seed_matches_rev = find_seed_matches(primer, rev_comp_chunk)
                for match_pos in seed_matches_rev:
                    if match_pos + len(primer) <= len(rev_comp_chunk):
                        full_match = rev_comp_chunk[match_pos:match_pos + len(primer)]
                        if get_hamming(primer, full_match) <= (len(primer) - threshold):
                            count += 1
                            total_count += 1
                            if total_count > max_occurrences:
                                break

                if total_count > max_occurrences:
                    break

            if total_count <= max_occurrences:
                filtered_primers.append((primer, start_pos))

        end_time = time.time()
        print(f"Primer specificity analysis took {end_time - start_time:.2f} seconds.")

        return filtered_primers


    def filter_primer_pairs(up_primers, down_primers, max_pairing=12):
        filtered_pairs = []

        for up_primer, up_start in up_primers:
            for down_primer, down_end in down_primers:
                pairing_count = len(up_primer) - get_hamming(up_primer, down_primer)
                if pairing_count <= max_pairing:
                    distance = down_end - up_start + 1
                    filtered_pairs.append((up_primer, down_primer, distance))

        return filtered_pairs


    def sanitize_filename(filename):

        return re.sub(r'[<>:"/\\|?*]', '_', filename)


    def save_to_excel_split(df, writer, sheet_name):
        max_rows = 1048576  # Excel row limit
        num_splits = (len(df) // max_rows) + 1
        for split_idx in range(num_splits):
            start_row = split_idx * max_rows
            end_row = start_row + max_rows
            df_split = df.iloc[start_row:end_row]
            df_split.to_excel(writer, sheet_name=f'{sheet_name}_{split_idx + 1}', index=False)


    def save_results_to_excel(results, seq1, seq2, seq5, seq6, mutation_desc, base_filename='crRNA_result'):
        if not results:
            print("No crRNA results found.")
            return


        processed_results = set()

        for i, result in enumerate(results, 1):

            unique_id = (result.get('PAM'), result.get('position_text'), result.get('protospacer'))
            if unique_id in processed_results:
                continue
            processed_results.add(unique_id)

            df = pd.DataFrame([result])
            pam = result.get('PAM', 'unknown')
            position = result.get('protospacer', i)
            output_file = f'{base_filename}_PAM {pam}_Type of mutation {mutation_desc}.xlsx'


            output_file = sanitize_filename(output_file)

            try:
                detailed_results, summary_results = perform_specificity_analysis(result['protospacer'], seq5)
                detailed_results_rev, summary_results_rev = perform_specificity_analysis(result['protospacer'], seq6)


                rna_structure_data = predict_rna_structure_for_guide_sequences(
                    result.get('guide sequence') or result.get('protospacer'))

                site = result['pos in sequence'] - 1
                up_primers, down_primers = generate_primers(seq1, site)
                filtered_up_primers = filter_primers(up_primers, seq5, seq6, result['protospacer'])
                filtered_down_primers = filter_primers(down_primers, seq5, seq6, result['protospacer'])
                filtered_pairs = filter_primer_pairs(filtered_up_primers, filtered_down_primers)
                primer_df = pd.DataFrame(filtered_pairs, columns=['Up-Primer', 'Down-Primer', 'Distance'])

                with pd.ExcelWriter(output_file) as writer:
                    df.to_excel(writer, sheet_name='crRNA Results', index=False)


                    filtered_detailed_results = [res for res in detailed_results if
                                                 res['Total Match Count'] >= 10 and res['First Five Match Count'] >= 4]
                    if filtered_detailed_results:
                        save_to_excel_split(pd.DataFrame(filtered_detailed_results), writer, 'Detailed Results Forward')

                    if summary_results:
                        save_to_excel_split(pd.DataFrame(summary_results), writer, 'Summary Forward')


                    filtered_detailed_results_rev = [res for res in detailed_results_rev if
                                                     res['Total Match Count'] >= 10 and res[
                                                         'First Five Match Count'] >= 4]
                    if filtered_detailed_results_rev:
                        save_to_excel_split(pd.DataFrame(filtered_detailed_results_rev), writer,
                                            'Detailed Results Reverse')

                    if summary_results_rev:
                        save_to_excel_split(pd.DataFrame(summary_results_rev), writer, 'Summary Reverse')

                    # Save RNA structure prediction
                    pd.DataFrame(rna_structure_data).to_excel(writer, sheet_name='RNA Structure', index=False)

                    # Save primers
                    save_to_excel_split(primer_df, writer, 'Primers')

                print(f"Result saved to {os.path.abspath(output_file)}")
            except Exception as e:
                print(f"Failed to save result {i}: {e}")


    def process_sequence(text):

        match_between = re.search(r"5'(.*?)3'", text, re.DOTALL)
        if match_between:
            sequence_between = re.sub(r'[^ATCG]', '', match_between.group(1))
        else:
            sequence_between = ""


        match_after = re.search(r"3'(.*)", text, re.DOTALL)
        if match_after:
            sequence_after = re.sub(r'[^ATCG]', '', match_after.group(1))
        else:
            sequence_after = ""


        combined_sequence = sequence_between + "H" + sequence_after

        return combined_sequence


    rsnumber = input("Please enter the RS number：")
    sequence5 = read_fasta_sequence(input("Please enter multiple FASTA file names (separated by spaces): ").split())
    sequence6 = reverse_complement(sequence5)

    options = webdriver.ChromeOptions()
    # options.add_argument('--headless')
    driver = webdriver.Chrome(service=Service(ChromeDriverManager().install()), options=options)

    try:
        url = f"https://www.ncbi.nlm.nih.gov/snp/{rsnumber}#flanks"
        driver.get(url)


        time.sleep(5)
        print("The page loads completely")


        select = WebDriverWait(driver, 30).until(
            EC.presence_of_element_located((By.ID, "flank_length"))
        )
        select.click()
        print("The selection box has been clicked")


        option_200nt = WebDriverWait(driver, 30).until(
            EC.element_to_be_clickable((By.XPATH, "//option[@value='200']"))
        )
        option_200nt.click()
        print("The 200 nt option is selected")


        retrieve_button = WebDriverWait(driver, 30).until(
            EC.element_to_be_clickable((By.ID, "retrieve_flank"))
        )
        retrieve_button.click()
        print("The search button has been clicked")


        print("Wait for the results to load...")
        time.sleep(10)


        sequence_text = driver.execute_script("return document.getElementById('flanking_sequence').innerText;")
        print("200 nt Flanking Sequence:", sequence_text)

        position_dt = WebDriverWait(driver, 30).until(
            EC.presence_of_element_located((By.XPATH, "//dt[text()='Position']"))
        )
        position_dd = position_dt.find_element(By.XPATH, "following-sibling::dd/span[1]")
        position_text = position_dd.text
        print("Position Information:", position_text)

    finally:
        driver.quit()


    text = sequence_text

    sequence1 = process_sequence(text)
    print("Processed sequences:", sequence1)

    try:
        rs_id = rsnumber.strip()
        alleles1, alleles2 = fetch_snp_data(rs_id)
        if alleles1:
            print(f"Alleles data: {alleles1}")


            replacement_char = alleles2[0]
            seq1 = sequence1[:200] + replacement_char + sequence1[201:]


            mutated_sequences = introduce_mutations(seq1, alleles1)

            all_results = []
            for seq_name, mutated_seq in mutated_sequences.items():

                results_forward, mutation_desc_forward = process_sequences(seq1, mutated_seq, "正向")


                sequence3 = reverse_complement(seq1)
                sequence4 = reverse_complement(mutated_seq)
                results_reverse, mutation_desc_reverse = process_sequences(sequence3, sequence4, "反向")

                all_results = results_forward + results_reverse


                save_results_to_excel(all_results, seq1, sequence3, sequence5, sequence6,
                                      mutation_desc_forward or mutation_desc_reverse)

            input("Press Enter to exit...")
        else:
            print("Unable to get the data, please check if the RS number is correct")

    except Exception as e:
        print(f"An error occurred: {e}")
        input("Press Enter to exit...")

if F == 'C':
    import time
    import re
    import pandas as pd
    import os
    import RNA
    import csv
    from Bio import SeqIO


    def fetch_snp_data(rs_id):
        alleles1 = []
        alleles2_merged = None


        current_directory = os.getcwd()


        csv_file_path = os.path.join(current_directory, 'GCF_000001405.40_1_5')


        try:
            with open(csv_file_path, 'r', encoding='utf-8') as file:
                reader = csv.reader(file, delimiter='\t')
                for line_num, line in enumerate(reader, start=1):

                    if not line:
                        continue


                    if len(line) < 5:
                        continue


                    if line[2] == rs_id:

                        alleles1 = line[4].split(',')
                        alleles2_merged = line[3]
                        break
            if alleles2_merged is None:
                print(f"SNP {rs_id} 未在文件中找到。")
                return None, None

            return alleles1, alleles2_merged

        except Exception as e:
            print(f"读取文件时发生错误：{e}")
            return None, None


    def introduce_mutations(seq, alleles, position=201):
        seqs_with_mutations = {}
        for i, allele in enumerate(alleles, 1):
            mutated_seq = list(seq)
            mutated_seq[position - 1] = allele
            seqs_with_mutations[f'seq2.{i}'] = ''.join(mutated_seq)

        return seqs_with_mutations


    def reverse_complement(sequence):
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        return ''.join(complement[base] for base in reversed(sequence.upper()))


    def find_differences(str1, str2):
        min_len = min(len(str1), len(str2))
        differences = []
        for i in range(min_len):
            if str1[i].lower() != str2[i].lower():
                differences.append(i)
        if len(str1) != len(str2):
            differences.extend(range(min_len, max(len(str1), len(str2))))
        return differences


    def read_fasta_sequence(file_names):
        sequence = []
        for file_name in file_names:
            with open(file_name, 'r') as file:
                for line in file:
                    if not line.startswith('>'):

                        cleaned_line = re.sub(r'[^ATCGN]', '', line.strip())
                        sequence.append(cleaned_line)
        return ''.join(sequence)


    def get_mutation_description(seq1, seq2, differences):
        mutations = []
        for pos in differences:
            mutations.append(f"{seq1[pos]}->{seq2[pos]} at {position_text}")
        return ', '.join(mutations)


    def mutate_string(g, pos):
        mutation_dict = {
            'A': ['C', 'T', 'G'],
            'T': ['C', 'A', 'G'],
            'G': ['C', 'A', 'T'],
            'C': ['A', 'G', 'T']
        }
        original_case = [(ch.isupper(), i) for i, ch in enumerate(g)]
        g_upper = g.upper()
        pos -= 1
        if len(g_upper) < 5:
            return [g]
        result = []
        for i in range(2, 5):
            if i == pos:
                continue
            char = g_upper[i]
            if char in mutation_dict:
                mutations = mutation_dict[char]
            else:
                return [g]
            for mutation in mutations:
                mutated_g = g_upper[:i] + mutation + g_upper[i + 1:]
                mutated_g_with_case = ''.join(
                    (ch.lower() if not ig_upper else ch.upper())
                    for ch, (ig_upper, idx) in zip(mutated_g, original_case)
                )
                result.append(mutated_g_with_case)
        return result


    def find_pattern_and_report(seq1, seq2, pos, direction, mutation_desc):
        results = []
        seq1 = seq1.upper()
        seq2 = seq2.upper()
        stand_symbol = '+' if direction == '正向结果' else '-'
        if seq1[pos] != 'T' and seq2[pos] == 'T':
            if pos < len(seq1) - 3 and seq1[pos + 1:pos + 3] in ["TA", "TC", "TG", "AT", "CT", "GT"]:
                pam = seq2[pos:pos + 4]
                crRNA_sequence = seq2[pos + 4:pos + 22]
                results.append({
                    'stand': stand_symbol,
                    'pos in sequence': pos + 1,
                    'PAM': pam,
                    'protospacer': crRNA_sequence,
                    'pos in PAM': 1,
                    'Description of the mutation': mutation_desc
                })
            if pos > 1 and seq1[pos - 2:pos] in ["TA", "TC", "TG", "AT", "CT", "GT"]:
                pam = seq2[pos - 2:pos + 2]
                crRNA_sequence = seq2[pos + 2:pos + 20]
                results.append({
                    'stand': stand_symbol,
                    'pos in sequence': pos + 1,
                    'PAM': pam,
                    'protospacer': crRNA_sequence,
                    'pos in PAM': 3,
                    'Description of the mutation': mutation_desc
                })
            if pos > 0 and pos < len(seq1) - 2 and seq1[pos - 1] == 'T' and seq1[pos + 1] in ["A", "C", "G"]:
                pam = seq2[pos - 1:pos + 3]
                crRNA_sequence = seq2[pos + 3:pos + 21]
                results.append({
                    'stand': stand_symbol,
                    'pos in sequence': pos + 1,
                    'PAM': pam,
                    'protospacer': crRNA_sequence,
                    'pos in PAM': 2,
                    'Description of the mutation': mutation_desc
                })
            if pos < len(seq1) - 1 and seq1[pos + 1] == 'T' and seq1[pos - 1] in ["A", "C", "G"]:
                pam = seq2[pos - 1:pos + 3]
                crRNA_sequence = seq2[pos + 3:pos + 21]
                results.append({
                    'stand': stand_symbol,
                    'pos in sequence': pos + 1,
                    'PAM': pam,
                    'protospacer': crRNA_sequence,
                    'pos in PAM': 2,
                    'Description of the mutation': mutation_desc
                })

        patterns = [r'ttt[atcg]', r't[atcg]t[atcg]', r'[atcg]tt[atcg]', r'tt[atcg][atcg]']
        start = max(0, pos - 14)
        context = seq2[start:pos]

        for pattern in patterns:
            i = 0
            while i <= len(context) - 4:
                match = re.match(pattern, context[i:], re.IGNORECASE)
                if match:
                    end_pos = match.end() + i + start
                    following_chars = seq2[end_pos:end_pos + 18]
                    new_pos = pos - end_pos + 1
                    result = mutate_string(following_chars, new_pos)
                    results.append({
                        'stand': stand_symbol,
                        'pos in sequence': pos + 1,
                        'PAM': match.group().upper(),
                        'protospacer': following_chars,
                        'pos in protospacer': new_pos,
                        'guide sequence': result,
                        'Description of the mutation': mutation_desc
                    })
                    i += 1
                else:
                    i += 1

        return results


    def process_sequences(seq1, seq2, direction):
        differences = find_differences(seq1, seq2)
        mutation_desc = get_mutation_description(seq1, seq2, differences)
        results = []
        for pos in differences:
            results.extend(find_pattern_and_report(seq1, seq2, pos, direction, mutation_desc))
        return results, mutation_desc


    def find_pam_sites(seq1):
        pam_patterns = ["TTTX", "XTTX", "TXTX", "TTXX"]
        pam_sites = []

        seq1 = seq1.upper()

        for i in range(len(seq1) - 3):
            pam_candidate = seq1[i:i + 4]
            for pattern in pam_patterns:
                if all(p == 'X' or p == c for p, c in zip(pattern, pam_candidate)):
                    pam_sites.append(i)
                    break

        return pam_sites


    def compare_crrna_with_targets(crrna, targets):
        results = []

        crrna = crrna.upper()

        for target in targets:
            target = target.upper()
            match_count = sum(1 for a, b in zip(crrna, target) if a == b)
            match_count_first_five = sum(1 for a, b in zip(crrna[:5], target[:5]) if a == b)
            results.append((match_count, match_count_first_five))
        return results


    def perform_specificity_analysis(crrna, seq1, chunk_size=100000, overlap=30):
        detailed_results = []
        summary_results = {}

        start_time = time.time()


        for start in range(0, len(seq1), chunk_size):
            end = min(start + chunk_size + overlap, len(seq1))
            chunk = seq1[start:end]


            pam_sites = find_pam_sites(chunk)
            targets = [(chunk[i:i + 4], chunk[i + 4:i + 22]) for i in pam_sites if i + 22 <= len(chunk)]

            for pam, target in targets:

                match_count_first_five = sum(1 for a, b in zip(crrna[:5], target[:5]) if a == b)

                if match_count_first_five < 4:
                    continue


                match_count_rest = sum(1 for a, b in zip(crrna[5:], target[5:]) if a == b)


                match_count = match_count_first_five + match_count_rest

                detailed_results.append({
                    'Index': len(detailed_results) + 1,
                    'PAM': pam,
                    'Target Sequence': target,
                    'Total Match Count': match_count,
                    'First Five Match Count': match_count_first_five
                })

                key = (match_count_first_five, match_count)
                if key not in summary_results:
                    summary_results[key] = 0
                summary_results[key] += 1

        summary_results = sorted(
            [{'First Five Match Count': k[0], 'Total Match Count': k[1], 'Count': v} for k, v in
             summary_results.items()],
            key=lambda x: (-x['First Five Match Count'], -x['Total Match Count'])
        )

        end_time = time.time()
        print(f"RNA specificity analysis took {end_time - start_time:.2f} seconds.")

        return detailed_results, summary_results


    def predict_rna_structure(sequence):
        fc = RNA.fold_compound(sequence)
        structure, mfe = fc.mfe()
        return structure, mfe


    def predict_rna_structure_for_guide_sequences(guide_sequences):
        rna_structure_data = []
        prefix_sequence = "TAATTTCTACTAAGTGTAGAT".replace("T", "U")


        if isinstance(guide_sequences, str):
            guide_sequences = [guide_sequences]


        for rna_sequence in guide_sequences:
            rna_sequence = rna_sequence.replace("T", "U")
            full_sequence = prefix_sequence + rna_sequence
            structure, mfe = predict_rna_structure(full_sequence)

            rna_structure_data.append({
                'Full RNA Sequence': full_sequence,
                'Predicted Structure': structure,
                'MFE (kcal/mol)': mfe
            })

        return rna_structure_data


    def get_hamming(seq1, seq2):
        if len(seq1) != len(seq2):
            raise ValueError("Sequences must be of equal length")
        return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


    def generate_primers(seq1, site, length=30):
        up_primers = []
        down_primers = []

        for i in range(max(0, site - 125), max(0, site - 50) - length + 1):
            up_primers.append((seq1[i:i + length].upper(), i))

        for i in range(site + 50, site + 125 - length):
            down_seq = seq1[i:i + length]
            down_primers.append((reverse_complement(down_seq).upper(), i + length - 1))

        return up_primers, down_primers


    def calculate_gc_content(seq):
        gc_count = seq.upper().count('G') + seq.upper().count('C')
        return gc_count / len(seq)


    def calculate_rna_pairing(primer, rna_seq):
        max_pairing = 0
        for i in range(len(rna_seq) - len(primer) + 1):
            pairing = len(primer) - get_hamming(primer, rna_seq[i:i + len(primer)])
            if pairing > max_pairing:
                max_pairing = pairing
        return max_pairing


    def filter_primers(primers, seq1, rev_comp_seq1, rna_seq, threshold=20, max_occurrences=22, gc_min=0.3, gc_max=0.7,
                       max_self_pairing=6, max_rna_pairing=6, seed_length=15, chunk_size=100000, overlap=30):
        filtered_primers = []

        start_time = time.time()

        def find_seed_matches(primer, sequence):
            seed = primer[:seed_length]
            positions = []
            for i in range(len(sequence) - seed_length + 1):
                if sequence[i:i + seed_length] == seed:
                    positions.append(i)
            return positions

        def contains_n(sequence):

            return 'N' in sequence

        for primer, start_pos in primers:
            gc_content = calculate_gc_content(primer)
            if gc_content < gc_min or gc_content > gc_max:
                continue

            reverse_comp = reverse_complement(primer)
            self_pairing_count = len(primer) - get_hamming(primer, reverse_comp)
            if self_pairing_count > max_self_pairing:
                continue

            rna_pairing_count = calculate_rna_pairing(primer, rna_seq)
            if rna_pairing_count > max_rna_pairing:
                continue

            total_count = 0

            for start in range(0, len(seq1), chunk_size):
                end = min(start + chunk_size + overlap, len(seq1))
                chunk = seq1[start:end]
                rev_comp_chunk = rev_comp_seq1[start:end]


                if contains_n(chunk) or contains_n(rev_comp_chunk):
                    continue

                count = 0
                seed_matches = find_seed_matches(primer, chunk)
                for match_pos in seed_matches:
                    if match_pos + len(primer) <= len(chunk):
                        full_match = chunk[match_pos:match_pos + len(primer)]
                        if get_hamming(primer, full_match) <= (len(primer) - threshold):
                            count += 1
                            total_count += 1
                            if total_count > max_occurrences:
                                break

                if total_count > max_occurrences:
                    break

                seed_matches_rev = find_seed_matches(primer, rev_comp_chunk)
                for match_pos in seed_matches_rev:
                    if match_pos + len(primer) <= len(rev_comp_chunk):
                        full_match = rev_comp_chunk[match_pos:match_pos + len(primer)]
                        if get_hamming(primer, full_match) <= (len(primer) - threshold):
                            count += 1
                            total_count += 1
                            if total_count > max_occurrences:
                                break

                if total_count > max_occurrences:
                    break

            if total_count <= max_occurrences:
                filtered_primers.append((primer, start_pos))

        end_time = time.time()
        print(f"Primer specificity analysis took {end_time - start_time:.2f} seconds.")

        return filtered_primers


    def filter_primer_pairs(up_primers, down_primers, max_pairing=12):
        filtered_pairs = []

        for up_primer, up_start in up_primers:
            for down_primer, down_end in down_primers:
                pairing_count = len(up_primer) - get_hamming(up_primer, down_primer)
                if pairing_count <= max_pairing:
                    distance = down_end - up_start + 1
                    filtered_pairs.append((up_primer, down_primer, distance))

        return filtered_pairs


    def sanitize_filename(filename):

        return re.sub(r'[<>:"/\\|?*]', '_', filename)


    def save_to_excel_split(df, writer, sheet_name):
        max_rows = 1048576  # Excel row limit
        num_splits = (len(df) // max_rows) + 1
        for split_idx in range(num_splits):
            start_row = split_idx * max_rows
            end_row = start_row + max_rows
            df_split = df.iloc[start_row:end_row]
            df_split.to_excel(writer, sheet_name=f'{sheet_name}_{split_idx + 1}', index=False)


    def save_results_to_excel(results, seq1, seq2, seq5, seq6, mutation_desc, base_filename='crRNA_result'):
        if not results:
            print("No crRNA results found.")
            return


        processed_results = set()

        for i, result in enumerate(results, 1):

            unique_id = (result.get('PAM'), result.get('position_text'), result.get('protospacer'))
            if unique_id in processed_results:
                continue
            processed_results.add(unique_id)

            df = pd.DataFrame([result])
            pam = result.get('PAM', 'unknown')
            position = result.get('protospacer', i)
            output_file = f'{base_filename}_PAM {pam}_Type of mutation {mutation_desc}.xlsx'


            output_file = sanitize_filename(output_file)

            try:
                detailed_results, summary_results = perform_specificity_analysis(result['protospacer'], seq5)
                detailed_results_rev, summary_results_rev = perform_specificity_analysis(result['protospacer'], seq6)

                # Predict RNA structure for the crRNA
                rna_structure_data = predict_rna_structure_for_guide_sequences(
                    result.get('guide sequence') or result.get('protospacer'))
                # Generate primers
                site = result['pos in sequence'] - 1
                up_primers, down_primers = generate_primers(seq1, site)
                filtered_up_primers = filter_primers(up_primers, seq5, seq6, result['protospacer'])
                filtered_down_primers = filter_primers(down_primers, seq5, seq6, result['protospacer'])
                filtered_pairs = filter_primer_pairs(filtered_up_primers, filtered_down_primers)
                primer_df = pd.DataFrame(filtered_pairs, columns=['Up-Primer', 'Down-Primer', 'Distance'])

                with pd.ExcelWriter(output_file) as writer:
                    df.to_excel(writer, sheet_name='crRNA Results', index=False)


                    filtered_detailed_results = [res for res in detailed_results if
                                                 res['Total Match Count'] >= 10 and res['First Five Match Count'] >= 4]
                    if filtered_detailed_results:
                        save_to_excel_split(pd.DataFrame(filtered_detailed_results), writer, 'Detailed Results Forward')

                    if summary_results:
                        save_to_excel_split(pd.DataFrame(summary_results), writer, 'Summary Forward')


                    filtered_detailed_results_rev = [res for res in detailed_results_rev if
                                                     res['Total Match Count'] >= 10 and res[
                                                         'First Five Match Count'] >= 4]
                    if filtered_detailed_results_rev:
                        save_to_excel_split(pd.DataFrame(filtered_detailed_results_rev), writer,
                                            'Detailed Results Reverse')

                    if summary_results_rev:
                        save_to_excel_split(pd.DataFrame(summary_results_rev), writer, 'Summary Reverse')

                    # Save RNA structure prediction
                    pd.DataFrame(rna_structure_data).to_excel(writer, sheet_name='RNA Structure', index=False)

                    # Save primers
                    save_to_excel_split(primer_df, writer, 'Primers')

                print(f"Result saved to {os.path.abspath(output_file)}")
            except Exception as e:
                print(f"Failed to save result {i}: {e}")


    def process_sequence(sequence_text):

        return sequence_text.strip()


    def get_sequence_from_fasta(genome_id, position, fasta_file_path):

        with open(fasta_file_path, 'r') as fasta_file:
            for record in SeqIO.parse(fasta_file, "fasta"):

                if genome_id in record.id:

                    sequence = str(record.seq)
                    start, end = position
                    if start >= 200 and end + 200 <= len(sequence):

                        left_seq = sequence[start - 200:start]
                        right_seq = sequence[end:end + 200]
                        return left_seq + 'H' + right_seq
                    else:
                        raise ValueError("This location is at the edge of the genome and Smart-SNPer:A is recommended")
        return None


    rsnumber = input("Please enter the RS number：")
    file_names_input = input(
        "Please enter multiple FASTA file names (separated by spaces), or press Enter to use default genomes : ").strip()


    if not file_names_input:
        file_names_input = "chrY.fasta chrX.fasta chrMT.fasta chr1.fasta chr2.fasta chr3.fasta chr4.fasta chr5.fasta chr6.fasta chr7.fasta chr8.fasta chr9.fasta chr10.fasta chr11.fasta chr12.fasta chr13.fasta chr14.fasta chr15.fasta chr16.fasta chr17.fasta chr18.fasta chr19.fasta chr20.fasta chr21.fasta chr22.fasta"


    file_names_list = file_names_input.split()
    sequence5 = read_fasta_sequence(file_names_list)
    sequence6 = reverse_complement(sequence5)

    fasta_file_path = "sequence.fasta"
    csv_file_path = "GCF_000001405.40_1_5"


    with open(csv_file_path, 'r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:

            if len(row) < 3:
                continue

            try:
                rsid = row[2]

                if rsnumber == rsid:
                    genome_id = row[0]
                    position = int(row[1])


                    try:
                        sequence1 = get_sequence_from_fasta(genome_id, (position, position + 1), fasta_file_path)
                        print("Processed sequences:", sequence1)
                        position_text = f"{genome_id}:{position}"
                        print("Position Information:", position_text)
                    except ValueError as e:
                        print(str(e))
                    break
            except ValueError as e:
                continue

    try:
        rs_id = rsnumber.strip()
        alleles1, alleles2 = fetch_snp_data(rs_id)
        if alleles1:
            print(f"Alleles data: {alleles1}")


            replacement_char = alleles2[0]
            seq1 = sequence1[:200] + replacement_char + sequence1[201:]


            mutated_sequences = introduce_mutations(seq1, alleles1)

            all_results = []
            for seq_name, mutated_seq in mutated_sequences.items():

                results_forward, mutation_desc_forward = process_sequences(seq1, mutated_seq, "正向")


                sequence3 = reverse_complement(seq1)
                sequence4 = reverse_complement(mutated_seq)
                results_reverse, mutation_desc_reverse = process_sequences(sequence3, sequence4, "反向")

                all_results = results_forward + results_reverse


                save_results_to_excel(all_results, seq1, sequence3, sequence5, sequence6,
                                      mutation_desc_forward or mutation_desc_reverse)

            input("Press Enter to exit...")
        else:
            print("Unable to get the data, please check if the RS number is correct")


    except Exception as e:
        print(f"An error occurred: {e}")
        input("Press Enter to exit...")
else:
    print('Typing error')
    input("Press Enter to exit...")