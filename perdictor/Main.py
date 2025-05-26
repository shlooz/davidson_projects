from SecondaryStructurePredictor import StructurePredictor

class Main:
    def __init__(self, seq):
        self.propensities_path = r"C:\Users\eitha\programing\Waitzman\final_project\propensities.txt"
        self.sequence = seq  # Example sequence
        self.window_size = 6
        self.min_region_length = 4
        self.min_propensity = 0.4  # Adjusted for normalized propensities

    def run(self):
        predictor = StructurePredictor(self.propensities_path)
        regions = predictor.predict(
            sequence=self.sequence,
            window_size=self.window_size,
            min_region_length=self.min_region_length,
            min_propensity=self.min_propensity
        )
        for start, end, struct in regions:
            print(f"Region: {start}-{end}, Structure: {struct}")

if __name__ == "__main__":
    code = Main("ARNDCEQGHILKMFPSTWYV")
    code.run()
