import uproot
import matplotlib.pyplot as plt
import mplhep as hep

hep.style.use("ROOT")
plt.style.use(hep.style.ROOT)

def plot(hName, xLabel, xMin, xMax):
    fig = plt.figure()
    ax = fig.subplots()

    for procName in procs:
        hist_tot = None
        for proc in procs[procName]:
            file = uproot.open(f"output/h_ww/histmaker/{proc}.root")
            hist = file[hName].to_hist()
            if hist_tot == None:
                hist_tot = hist
            else:
                hist_tot += hist
        hist_tot /= sum(hist_tot.values())
        hep.histplot(hist_tot, label=procName, ax=ax)

    ax.legend(fontsize='x-small')
    ax.set_xlabel(xLabel)
    ax.set_ylabel("Events")
    ax.set_xlim(xMin, xMax)
    #ax.set_yscale('log')

    plt.savefig(f"{outDir}/{hName}.png", bbox_inches="tight")
    plt.close()



if __name__ == "__main__":
    outDir = "/home/submit/jaeyserm/public_html/fccee/h_ww/"

    procs = {}
    procs["qqHWW"] = ['wzp6_ee_ssH_HWW_ecm240', 'wzp6_ee_ccH_HWW_ecm240', 'wzp6_ee_bbH_HWW_ecm240', 'wzp6_ee_qqH_HWW_ecm240']
    procs["nunuHWW"] = ['wzp6_ee_nunuH_HWW_ecm240']
    procs["eeHWW"] = ['wzp6_ee_eeH_HWW_ecm240']
    procs["mumuHWW"] = ['wzp6_ee_mumuH_HWW_ecm240']
    procs["tautauHWW"] = ['wzp6_ee_tautauH_HWW_ecm240']

    plot("missingMass", "Missing mass (GeV)", 0, 250)
    plot("electrons_all_p", "Electrons momentum (GeV)", 0, 100)
    plot("muons_all_p", "Muons momentum (GeV)", 0, 100)




