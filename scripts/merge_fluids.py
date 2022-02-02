from netprop.models import PropagationResultModel


def merge_all_liquids(res: PropagationResultModel, unified_liquid_name: str="info"):
    for n in res.nodes:
        res[n].liquids[unified_liquid_name] = sum(res[n].liquids.values(), 0)
