import uproot
import matplotlib.pyplot as plt
import mplhep as hep

hep.style.use("ROOT")
plt.style.use(hep.style.ROOT)

def plot(fIn, fOut):
    file = uproot.open(fIn)
    score_B = file["recojet_isB"].to_hist()
    score_C = file["recojet_isC"].to_hist()
    score_S = file["recojet_isS"].to_hist()

    fig = plt.figure()
    ax = fig.subplots()

    hep.histplot(score_B, label="b-score", ax=ax)
    hep.histplot(score_C, label="c-score", ax=ax)
    hep.histplot(score_S, label="s-score", ax=ax)

    ax.legend(fontsize='x-small')
    ax.set_xlabel("Tagger score")
    ax.set_ylabel("Events")
    ax.set_xlim(0, 1)
    ax.set_yscale('log')

    plt.savefig(f"{fOut}.png", bbox_inches="tight")
    plt.close()  # Close the figure to avoid display in some environments



if __name__ == "__main__":
    plot("output/clustering_tagging/histmaker/wzp6_ee_nunuH_Hbb_ecm240.root", "nunuH_Hbb")
    plot("output/clustering_tagging/histmaker/wzp6_ee_nunuH_Hcc_ecm240.root", "nunuH_Hcc")
    plot("output/clustering_tagging/histmaker/wzp6_ee_nunuH_Hss_ecm240.root", "nunuH_Hss")