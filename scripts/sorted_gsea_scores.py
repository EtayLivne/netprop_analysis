from utils.utils import load_json

scores = load_json(r"C:\studies\gsea_scores.json")

for k in sorted(list(scores.keys()), key=lambda k: len(scores[k][1]), reverse=False)[:10]:
    print(f"{k}: {scores[k][0]}, {len(scores[k][1])}")