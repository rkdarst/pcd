# Richard Darst, July 2011

from math import log, exp

initialTime = 1
# 1.25892541179417 = exp(log(10) / 10)
logConstant = 1.25892541179417
logInterval = 10.
def logTime(timeIndex=1):
    return (initialTime * logConstant**timeIndex)
def logTimeIndex(time):
    return log(time/initialTime) / (log(10)/10)

log2 = lambda x: log(x, 2)

def mutual_information(G0, G1):
    assert G0.n == G1.n
    N = G0.n

    MI = 0.0
    for c0 in [    c for c in range(G0.Ncmty) if G0.cmtyN[c] != 0]:
        for c1 in [c for c in range(G1.Ncmty) if G1.cmtyN[c] != 0]:
            # We correlate c0 in G0 and c1 in G1 according to the formula.
            n0 = G0.cmtyN[c0]
            n1 = G1.cmtyN[c1]

            # number of shared particles?
            n_shared = 0
            for n in    G0.cmtyll[c0, :n0]:
                if n in G1.cmtyll[c1, :n1]:
                    n_shared += 1
            s0 = set(G0.cmtyll[c0, :n0])
            s1 = set(G1.cmtyll[c1, :n1])
            assert n_shared == len(s0 & s1)

            if n_shared == 0:
                continue

            MI += (n_shared/float(N)) * log2(n_shared*N/float(n0*n1))

    return MI

if __name__ == "__main__":
    print logTime(-10)
    print logTime(logTimeIndex(100))
