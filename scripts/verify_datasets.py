#!/usr/bin/env python3
"""01-pre: verify candidate GEO accessions via NCBI E-utilities (db=gds).
Pulls n_samples, platform(s), organism, type, title for each accession so we
can decide HF main/validation cohorts and whether the sarcopenia side can
support WGCNA (rule of thumb: per-group N > 15-20)."""
import requests, time, sys, json

EMAIL = "candicewu0515@gmail.com"
BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

CANDIDATES = {
    "HF / LV (heart)": ["GSE57338", "GSE5406", "GSE16499", "GSE1145", "GSE76701", "GSE116250"],
    "Sarcopenia / aging muscle": ["GSE8479", "GSE1428", "GSE25941", "GSE9103", "GSE38718", "GSE111016", "GSE167186"],
    "RED-LINE (avoid)": ["GSE56815"],
}

def esearch_acc(acc):
    r = requests.get(f"{BASE}/esearch.fcgi", params={
        "db": "gds", "term": f"{acc}[ACCN] AND gse[ETYP]",
        "retmode": "json", "email": EMAIL}, timeout=30)
    r.raise_for_status()
    return r.json()["esearchresult"]["idlist"]

def esummary(uids):
    r = requests.post(f"{BASE}/esummary.fcgi", data={
        "db": "gds", "id": ",".join(uids), "retmode": "json", "email": EMAIL}, timeout=30)
    r.raise_for_status()
    return r.json()["result"]

def main():
    out = []
    for group, accs in CANDIDATES.items():
        print(f"\n{'='*78}\n{group}\n{'='*78}")
        for acc in accs:
            try:
                uids = esearch_acc(acc)
                time.sleep(0.4)
                if not uids:
                    print(f"  {acc:12s}  NOT FOUND")
                    out.append({"group": group, "acc": acc, "status": "NOT_FOUND"})
                    continue
                res = esummary(uids); time.sleep(0.4)
                # pick the exact GSE record
                rec = None
                for uid in res.get("uids", []):
                    s = res[uid]
                    if s.get("accession", "").upper() == acc.upper():
                        rec = s; break
                if rec is None and res.get("uids"):
                    rec = res[res["uids"][0]]
                n = rec.get("n_samples")
                plat = rec.get("gpl", "")
                taxon = rec.get("taxon", "")
                gdstype = rec.get("gdstype", "")
                title = (rec.get("title", "") or "")[:70]
                pdat = rec.get("pdat", "")
                print(f"  {acc:12s}  N={str(n):>4}  GPL={plat:<8}  {taxon:<13}  {pdat}")
                print(f"               {title}")
                print(f"               type: {gdstype}")
                out.append({"group": group, "acc": acc, "status": "OK", "n_samples": n,
                            "gpl": plat, "taxon": taxon, "gdstype": gdstype,
                            "title": rec.get("title", ""), "pdat": pdat,
                            "summary": (rec.get("summary", "") or "")[:400]})
            except Exception as e:
                print(f"  {acc:12s}  ERROR: {e}")
                out.append({"group": group, "acc": acc, "status": f"ERROR:{e}"})
    with open("results/00_config/dataset_verification.json", "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved -> results/00_config/dataset_verification.json ({len(out)} records)")

if __name__ == "__main__":
    main()
