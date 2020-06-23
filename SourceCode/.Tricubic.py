import numpy as np
def AdaptiveLPF(m, srt):
    derme = dfdt(m, srt)
    normVal = np.max(derme)
    badPointsExist = len(badPoints(derme, normVal, srt, len(m))) > 0

    while (badPointsExist):
        derme = dfdt(m, srt)
        a = badPoints(derme, normVal, srt, len(m))

        if len(a) == 0:
            badPointsExist = False
            break

        # print("--->")
        # print(len(a))
        m = dfdtFix(m, a, srt, 2)

    return m


def dfdt(m, srt):
    der2 = np.zeros(len(m))
    for i in range(1, len(m) - 1):
        der2[i] = (m[i + 1] - 2 * m[i] + m[i - 1])

    return der2


def badPoints(der2, val, srt, limit):
    der2 = der2 / val
    if srt == 0.1:
        tol = 0.03
    elif srt == 1:
        tol = 0.0005
    elif srt == 16:
        tol = 0.0005
    a = []
    for i in range(len(der2)):
        if der2[i] > tol and i > 50 and i < limit - 50:
            a.append(i)
    return a


def dfdtFix(m, a, srt, numPasses):
    # kernel=[1/3,1/3,1/3]
    kernel = [1 / 4, 2 / 4, 1 / 4]
    # kernel=[1/1048,11/1048,55/1048,165/1048,330/1048,462/1048,462/1048,330/1048,165/1048,55/1048,11/1048,1/2048]
    # kernel=[0.006,0.061,0.242,0.383,0.242,0.061,0.006]
    # kernel=[1/16,4/16,6/16,4/16,1/16]
    halfwidth = len(kernel) // 2
    window = 30

    for j in range(numPasses):
        for k in range(len(a)):
            fil = m
            for i in range(a[k] - window, a[k] + window):
                pixel = 0
                for p in range(len(kernel)):
                    pixel += m[i + p - halfwidth] * kernel[p]
                fil[i] = pixel
            m = fil
    return m