from pipelines.infection_resiliency_pipeline.run_pipeline import run_once


def main() -> None:
    individual_interactors_metadata = "/data/metadata.json"
    prop_results_file = "/data/all.csv"
    resiliency_dir = "/data/output"
    ratios = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    min_subset_size = 5
    max_subset_size = 30
    top_propagated_rank_threshold = 1000
    top_resilient_rank_threshold = 100
    run_name = "default"
    
    run_once(individual_interactors_metadata, prop_results_file, resiliency_dir,
             ratios, min_subset_size, max_subset_size, top_propagated_rank_threshold, top_resilient_rank_threshold,
             run_name, reset_state=False)




if __name__ == "__main__":
    main()