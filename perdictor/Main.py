import os
import statistics
import matplotlib.pyplot as plt
from SecondaryStructurePredictor import StructurePredictor

class Main:
    def __init__(self, pdb_path: str = ""):
        # Initialize parameters and default test folder
        self.pdb_path = pdb_path
        self.window_size = 4
        self.min_region_length = 2
        self.min_propensity = 0.1
        self.test_folder = r"C:\Users\eitha\programing\Waitzman\final_project\testFiles"

    def execute(self):
        # Run prediction and print results for a single PDB file
        predictor = StructurePredictor()
        predictor.run_prediction_from_pdb(
            pdb_path=self.pdb_path,
            window_size=self.window_size,
            min_region_length=self.min_region_length,
            min_propensity=self.min_propensity
        )

    def evaluate_all_files(self, isFoundPropensities: bool=True):
        # Evaluate predictions across all PDB files in test folder
        predictor = StructurePredictor()
        accuracies = []
        helix_f1s = []
        sheet_f1s = []

        for filename in os.listdir(self.test_folder):
            if filename.lower().endswith('.pdb'):
                full_path = os.path.join(self.test_folder, filename)
                print(f"\nEvaluating {filename}...")
                sequence = predictor.extract_sequence_from_pdb(full_path)
                if not sequence:
                    print(f"  Skipping {filename}: No sequence found.")
                    continue

                # Predict structures and get actual regions
                predicted = predictor.predict(
                    sequence,
                    window_size=self.window_size,
                    min_region_length=self.min_region_length,
                    min_propensity=self.min_propensity,
                    foundPropensities=isFoundPropensities
                )

                real = predictor.extract_real_structure_regions_from_pdb(full_path)

                # Score prediction
                acc = predictor.score_predictions(predicted, real, len(sequence))
                helix_f1 = predictor.score_helix_predictions(predicted, real, len(sequence))
                sheet_f1 = predictor.score_sheet_predictions(predicted, real, len(sequence))

                # Store scores
                accuracies.append(acc)
                helix_f1s.append(helix_f1)
                sheet_f1s.append(sheet_f1)

                print(f"  Accuracy: {acc * 100:.2f}%")
                print(f"  F1-score (Helix): {helix_f1:.2f}")
                print(f"  F1-score (Sheet): {sheet_f1:.2f}")

        if not accuracies:
            print("\nNo valid .pdb files with sequences found.")
            return

        # Compute statistics for accuracy and F1-scores
        avg_acc = statistics.mean(accuracies)
        med_acc = statistics.median(accuracies)
        max_acc = max(accuracies)
        min_acc = min(accuracies)

        avg_helix = statistics.mean(helix_f1s)
        med_helix = statistics.median(helix_f1s)
        max_helix = max(helix_f1s)
        min_helix = min(helix_f1s)

        avg_sheet = statistics.mean(sheet_f1s)
        med_sheet = statistics.median(sheet_f1s)
        max_sheet = max(sheet_f1s)
        min_sheet = min(sheet_f1s)

        # Print summary of results
        print(f"\n=== Summary for {len(accuracies)} PDB files ===")
        print(f"Average Accuracy:      {avg_acc * 100:.2f}%")
        print(f"Median Accuracy:       {med_acc * 100:.2f}%")
        print(f"Highest Accuracy:      {max_acc * 100:.2f}%")
        print(f"Lowest Accuracy:       {min_acc * 100:.2f}%")

        print(f"Average F1 (Helix):    {avg_helix * 100:.2f}")
        print(f"Median F1 (Helix):     {med_helix * 100:.2f}")
        print(f"Highest F1 (Helix):    {max_helix * 100:.2f}")
        print(f"Lowest F1 (Helix):     {min_helix * 100:.2f}")

        print(f"Average F1 (Sheet):    {avg_sheet * 100:.2f}")
        print(f"Median F1 (Sheet):     {med_sheet * 100:.2f}")
        print(f"Highest F1 (Sheet):    {max_sheet * 100:.2f}")
        print(f"Lowest F1 (Sheet):     {min_sheet * 100:.2f}")

        # Visualize results with a bar chart
        self.display_summary_chart(
            avg_acc, med_acc, max_acc, min_acc,
            avg_helix, med_helix, max_helix, min_helix,
            avg_sheet, med_sheet, max_sheet, min_sheet
        )

    def display_summary_chart(
        self,
        avg_acc, med_acc, max_acc, min_acc,
        avg_helix, med_helix, max_helix, min_helix,
        avg_sheet, med_sheet, max_sheet, min_sheet
    ):
        # Plot summary of accuracy and F1-scores as grouped bar chart
        categories = ["Total Accuracy", "Alpha Helix", "Beta Sheet"]
        averages = [avg_acc * 100, avg_helix * 100, avg_sheet * 100]
        medians = [med_acc * 100, med_helix * 100, med_sheet * 100]
        highs = [max_acc * 100, max_helix * 100, max_sheet * 100]
        lows = [min_acc * 100, min_helix * 100, min_sheet * 100]

        x = range(len(categories))
        width = 0.2

        plt.figure(figsize=(12, 6))
        plt.bar([i - 1.5 * width for i in x], averages, width=width, label='Average', color='cornflowerblue')
        plt.bar([i - 0.5 * width for i in x], medians, width=width, label='Median', color='seagreen')
        plt.bar([i + 0.5 * width for i in x], highs, width=width, label='Highest', color='gold')
        plt.bar([i + 1.5 * width for i in x], lows, width=width, label='Lowest', color='salmon')

        # Add text labels on bars
        for i, (avg, med, hi, lo) in enumerate(zip(averages, medians, highs, lows)):
            plt.text(i - 1.5 * width, avg + 1, f"{avg:.1f}", ha='center', fontsize=8)
            plt.text(i - 0.5 * width, med + 1, f"{med:.1f}", ha='center', fontsize=8)
            plt.text(i + 0.5 * width, hi + 1, f"{hi:.1f}", ha='center', fontsize=8)
            plt.text(i + 1.5 * width, lo + 1, f"{lo:.1f}", ha='center', fontsize=8)

        plt.ylabel("Score (%)")
        plt.title("Secondary Structure Prediction Summary")
        plt.xticks(ticks=x, labels=categories)
        plt.ylim(0, 110)
        plt.legend()
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        plt.tight_layout()
        plt.show()

    def showStructurePred(self, PDB_path, isFoundPropensities: bool=True):
        # Visualize structure prediction for a single PDB file
        predictor = StructurePredictor()
        structure = predictor.get_assignment_tuples(predictor.extract_sequence_from_pdb(PDB_path), isFoundPropensities)
        self.plot_amino_acid_structure(structure)

    def plot_amino_acid_structure(self, sequence):
        """
        Plot amino acid sequence colored by predicted structure.
        H=red, E=blue, C=black, gaps=gray.
        Adds residue index at start and end of each row.
        """

        color_map = {'H': 'red', 'E': 'blue', 'C': 'black', '-': 'gray'}



        # Create dictionary from position to (amino acid, structure)
        seq_dict = {pos: (aa, ss) for aa, ss, pos in sequence}

        min_pos = min(seq_dict.keys())
        max_pos = max(seq_dict.keys())

        # Fill full sequence including gaps
        full_sequence = []
        for pos in range(min_pos, max_pos + 1):
            if pos in seq_dict:
                aa, ss = seq_dict[pos]
            else:
                aa, ss = '-', '-'
            full_sequence.append((aa, ss, pos))

        # Plotting
        line_height = -1
        chars_per_line = 20
        total_lines = (len(full_sequence) + chars_per_line - 1) // chars_per_line

        plt.figure(figsize=(14, total_lines + 3))


        for i, (aa, ss, pos) in enumerate(full_sequence):
            line = i // chars_per_line
            x = (i % chars_per_line) + 2  # +2 to make room for index at start
            y = line_height * line
            color = color_map.get(ss, 'gray')
            plt.text(x, y, aa, fontsize=12, color=color, ha='center', va='center')

            # Write residue index at start of line
            if i % chars_per_line == 0:
                plt.text(0, y, str(pos), fontsize=10, color='dimgray', ha='right', va='center')

            # Write residue index at end of line
            if (i + 1) % chars_per_line == 0 or (i + 1) == len(full_sequence):
                plt.text(chars_per_line + 2.5, y, str(pos), fontsize=10, color='dimgray', ha='left', va='center')


                        # Colored title (H/E/C in color)
        plt.text(10, line_height * (total_lines + 1), "Amino Acid Sequence by Structure", fontsize=16, ha='center')
        plt.text(5.5, line_height * (total_lines + 2), "(", fontsize=14)
        plt.text(6.0, line_height * (total_lines + 2), "H", fontsize=14, color='red')
        plt.text(6.3, line_height * (total_lines + 2), "=red,", fontsize=14)
        plt.text(7.5, line_height * (total_lines + 2), "E", fontsize=14, color='blue')
        plt.text(7.8, line_height * (total_lines + 2), "=blue,", fontsize=14)
        plt.text(9.3, line_height * (total_lines + 2), "C", fontsize=14, color='black')
        plt.text(9.6, line_height * (total_lines + 2), "=black,", fontsize=14)
        plt.text(11.3, line_height * (total_lines + 2), "gaps=gray)", fontsize=14)

        # Format plot
        plt.xlim(-1, chars_per_line + 5)
        plt.ylim(line_height * (total_lines + 3), 1)
        plt.axis('off')
        plt.tight_layout()
        plt.show()
        


if __name__ == "__main__":
    pdb_path = r"C:\Users\eitha\programing\Waitzman\final_project\testFiles\1D4T.pdb"
    code = Main(pdb_path)
    code.execute()
    code.evaluate_all_files(True)
    code.showStructurePred(pdb_path, True)
