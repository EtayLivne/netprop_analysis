import pandas as pd
from functools import reduce
from pathlib import Path
import requests


def extract_index_sets(dfs: pd.DataFrame, threshold: int=100) -> list[set]:
    return [set(df.iloc[:threshold].index) for df in dfs]


def load_dfs(root_path: str) -> list[pd.DataFrame]:
    l = []
    for f in [str(f) for f in Path(root_path).glob("*")]:
        try:
            l.append(pd.read_csv(f, index_col="nodes"))
        except:
            print(f)
    return  l
    # return [pd.read_csv(str(f)) for f in Path(root_path).glob("*")]

#
# def intersect(index_sets: list[set]) -> set:
#     return reduce(lambda a, b: a & b, index_sets)


def top_intersected(index_sets: list[set]) -> dict:
    all_indexes = reduce(lambda a, b: a | b, index_sets)
    return {idx: len([_ for idx_set in index_sets if idx in idx_set]) for idx in all_indexes}


def normalize_and_threshold_score_dict(score_dict: dict[int, int], norm_factor: int, threshold: float=0.3):
    sorted_keys = sorted(list(score_dict.keys()), key=lambda k: score_dict[k], reverse=True)
    thresholded_keys = [k for k in sorted_keys if score_dict[k] >= threshold]
    #recreate dict
    d = {}
    for k in thresholded_keys:
        d[k] = score_dict[k] / norm_factor

    return d


def query_gestalt_ora(gene_list_file: str):
    webgestalt_url = "www.webgestalt.org/option.php"
    params = {
        "organism": "hsapiens",
        "enrich_method": "ORA"
    }


def analyze_by_mean_offset(top_ranked_threshold: int, result_dir_path: str, output_file_path: str) -> None:
    try:
        dfs = load_dfs(result_dir_path)
    except Exception as e:
        raise ValueError
    index_sets = extract_index_sets(dfs, threshold=top_ranked_threshold)

    # all = intersect(index_sets)
    ti = top_intersected(index_sets)
    normalized = normalize_and_threshold_score_dict(ti, len(index_sets))
    with open(output_file_path, "w") as handler:

        for k in normalized.keys():
            handler.write(str(k)+"\n")
    x = 7




# def main():
#     root_path = r"/pipelines/infection_resiliency_pipeline/pipeline_objs/temp_resiliency"
#     threshold = 100
#
#     try:
#         dfs = load_dfs(root_path)
#     except Exception as e:
#         raise ValueError
#     index_sets = extract_index_sets(dfs, threshold=threshold)
#
#     all = intersect(index_sets)
#     ti = top_intersected(index_sets)
#     normalized = normalize_and_threshold_score_dict(ti, len(index_sets))
#     with open("/gene_list.txt", "w") as handler:
#
#         for k in normalized.keys():
#             handler.write(str(k)+"\n")
#     x = 7
#
# if __name__ == "__main__":
#     # main()
#     # import pandas as pd
#     s = "http://www.webgestalt.org/process.php?organism=hsapiens&enrich_method=ORA&fdr_m \
#     ethod=BY&enriched_database_category=geneontology&enriched_database_name=Bio \
#     logical_Process_noRedundant&sig_method=top&sig_value=0.01&max_num=200&id_ty \
#     pe=entrezgene&gene_list=BMP2%0AAPC&ref_set=genome"
#     # res = requests.post(s)
#     # t = pd.read_html(res.content)
#     res = requests.post(s)
#     x = 7

