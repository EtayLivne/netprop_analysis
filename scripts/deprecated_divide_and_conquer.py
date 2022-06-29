# Returns all nodes in the graph that are not covid nodes themselves but have a covid neighbor
def get_covid_interactors(network: nx.Graph) -> set[str]:
    cov_interactors = set()
    for _, interactors in gen_covid_interactors_by_protein(network):
        cov_interactors |= interactors
    return cov_interactors


def gen_covid_interactors_by_protein(network: nx.Graph) -> set[str]:
    cov_nodes = [n for n, data in network.nodes(data=True) if data["species_id"] == "sars-cov-2"]
    for node in cov_nodes:
        yield node, set(network.neighbors(node))



def propagate_stefan():
    propagate_from_config(r"D:\configurations\stefan1\stefan_conf.json", ordering={"prior_set": 100})


def stefan_propagations(conf_files: list[str]):
    for file in conf_files:
        print(f"propagating {file}")
        propagate_from_config(str(file), ordering={"prior_set": 100}, max_processes=7   )


def intersect_top_propagated(results: list[str], k: int=50, output_file: str=None):
    sorted_results = sorted(results)    #alphabatical sorting would place the files as minor-major all the way
    majorities = sorted_results[::2]
    minorities = sorted_results[1::2]
    names = [res.split("\\")[-1].split("_majority")[0] for res in majorities]
    res_dict = dict()
    for i in range(len(majorities)):
        major_res = PropagationResultModel.parse_file(majorities[i])
        major_res_top_propagated = sorted([n for n in major_res.nodes.keys()], key=lambda n: major_res.nodes[n].liquids["info"], reverse=True)[:k]
        minor_res = PropagationResultModel.parse_file(minorities[i])
        minor_res_top_propagated = sorted([n for n in minor_res.nodes.keys()], key=lambda n: minor_res.nodes[n].liquids["info"], reverse=True)[:k]
        intersection = set(major_res_top_propagated) & set(minor_res_top_propagated)
        if len(intersection)/k > 0:
            print(f"{names[i]}: {len(intersection)/k}")
        res_dict[names[i]] = len(intersection)/k

    if output_file:
        with open(output_file, 'w') as handler:
            json.dump(res_dict, handler, indent=4)


def average_intersection(result_folders: list[str], k: int=50):
    scores = []
    for result_folder in result_folders:
        files = sorted([str(f) for f in Path(result_folder).glob("*.json")]) #alphabatical sorting would place the files as minor-major all the way
        majorities = files[::2]
        minorities = files[1::2]
        names = [res.split("\\")[-1].split("_majority")[0] for res in majorities]
        for i in range(len(majorities)):
            major_res = PropagationResultModel.parse_file(majorities[i])
            major_res_top_propagated = sorted([n for n in major_res.nodes.keys()], key=lambda n: major_res.nodes[n].liquids["info"], reverse=True)[:k]
            minor_res = PropagationResultModel.parse_file(minorities[i])
            minor_res_top_propagated = sorted([n for n in minor_res.nodes.keys()], key=lambda n: minor_res.nodes[n].liquids["info"], reverse=True)[:k]
            intersection = set(major_res_top_propagated) & set(minor_res_top_propagated)
            scores.append(len(intersection)/k)

    return scores


def top_prop_by_source(res_file: str, k: int=50):
    nodes = PropagationResultModel.parse_file(res_file).nodes
    liquids = nodes['1'].liquids.keys()
    return {
        liquid: sorted([n for n in nodes.keys()], key=lambda n: nodes[n].liquids[liquid], reverse=True)[:k]
        for liquid in liquids
    }


def most_intersected(dict_by_liquid: dict):
    nodes = {}
    for liquid, top_nodes in dict_by_liquid.items():
        for node in top_nodes:
            if node not in nodes:
                nodes[node] = []
            nodes[node].append(liquid)
    return nodes


def compare_results(res_1: str, res_2: str) -> bool:
    res_1_nodes = PropagationResultModel.parse_file(res_1).nodes
    res_2_nodes = PropagationResultModel.parse_file(res_2).nodes
    from math import fabs
    diffs = []
    counters = {"normal": 0, "weird": 0}
    for node, data in res_1_nodes.items():
        score_1, score_2 = data.liquids["info"], res_2_nodes[node].liquids["info"]

        if min(score_1, score_2) == 0:
            deviation = max(score_1, score_2)
            counters["weird"] += 1
        else:
            deviation = max(score_1, score_2) / min(score_1, score_2)
            counters["normal"] += 1

        #diffs.append(fabs(score_1- score_2)/max(score_1, score_2))
        diffs.append(deviation)
    print(f"max: {max(diffs)}\navg: {sum(diffs)/len(diffs)}")
    print(len([d for d in diffs if d < 1]))


# def validate_statistical_significance()


def more_detailed_intersection_data(results: list[str], k: int=50, output_file: str=None):
    sorted_results = sorted(results)    #alphabatical sorting would place the files as minor-major all the way
    majorities = sorted_results[::2]
    minorities = sorted_results[1::2]
    names = [res.split("\\")[-1].split("_majority")[0] for res in majorities]
    res_dict = {"chunks": dict(), "proteins": dict()}
    all_intersections = set()
    for i in range(len(majorities)):
        major_res = PropagationResultModel.parse_file(majorities[i])
        major_res_top_propagated = sorted([n for n in major_res.nodes.keys()], key=lambda n: major_res.nodes[n].liquids["info"], reverse=True)[:k]
        minor_res = PropagationResultModel.parse_file(minorities[i])
        minor_res_top_propagated = sorted([n for n in minor_res.nodes.keys()], key=lambda n: minor_res.nodes[n].liquids["info"], reverse=True)[:k]
        intersection = set(major_res_top_propagated) & set(minor_res_top_propagated)
        if len(intersection)/k > 0:
            print(f"{names[i]}: {len(intersection)/k}")
        res_dict["chunks"][names[i]] = {"size": len(intersection)/k, "intersecting_proteins": list(intersection)}
        all_intersections = all_intersections.union(intersection)

    res_dict["proteins"] = {
        p: len([1 for chunk in res_dict["chunks"].values() if p in chunk["intersecting_proteins"]]) for p in all_intersections
    }

    if output_file:
        with open(output_file, 'w') as handler:
            json.dump(res_dict, handler, indent=4)