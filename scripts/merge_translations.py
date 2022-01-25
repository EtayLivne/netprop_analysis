import json

def merge_translations():
    with open(r"D:\hgnc_complete.json", 'r') as handler:
        hgnc = json.load(handler)

    with open("D:\symbol_to_entrezgene_2021.json", 'r') as handler:
        my_translations = json.load(handler)

    for item in hgnc.values():
        symbol, prev_symbol, entrez_id = item["symbol"], item["prev_symbol"], item["entrez_id"]
        for i, sym in enumerate((symbol, prev_symbol)):
            if not sym:
                continue
            if sym not in my_translations:
                my_translations[sym] = entrez_id
            elif my_translations[sym] != entrez_id:
                print(f"hgnc sym {i}:{sym} clashes with my translation")
    with open(r"D:\complete_translation.json", 'w') as handler:
        json.dump(my_translations, handler, indent=4)


def decmopose_concatanated_symbols():
    with open(r"D:\complete_translation.json", 'r') as handler:
        hgnc = json.load(handler)
    new_dict = dict()
    for k, v in hgnc.items():
        if "|" in k:
            new_keys = k.split("|")
            for nk in new_keys:
                if nk in hgnc and hgnc[nk] != v:
                    pass
                else:
                    new_dict[nk] = v
    hgnc.update(new_dict)
    with open(r"D:\complete_translation_new.json", 'w') as handler:
        json.dump(hgnc, handler, indent=4)


if __name__ == "__main__":
    decmopose_concatanated_symbols()