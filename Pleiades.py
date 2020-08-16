# importing the necessary modules
import numpy as np
import matplotlib.pyplot as plt

# importing and storting the data for the Pleiades Cluster
cluster = 'Pleiades.csv'
cdata = np.loadtxt(cluster, delimiter=',', comments='#', usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16))

sourceid = cdata[:, 0]
ra = cdata[:, 1]
ra_err = cdata[:, 2]
dec = cdata[:, 3]
dec_err = cdata[:, 4]
para = cdata[:, 5]
para_err = cdata[:, 6]
pmra = cdata[:, 7]
pmra_err = cdata[:, 8]
pmdec = cdata[:, 9]
pmdec_err = cdata[:, 10]
gfluxovererr = cdata[:, 11]
gmeanmag = cdata[:, 12]
bprpexcess = cdata[:, 13]
bprp = cdata[:, 14]
agval = cdata[:, 15]
ebpminrpval = cdata[:, 16]

# Plotting the proper motion of all stars in the cluster data
plt.plot(pmra, pmdec, 'o', markersize = 0.1)
plt.xlim(-75, 100)
plt.ylim(-100, 75)
plt.xlabel('Proper Motion of Stars in Right Ascension')
plt.ylabel('Proper Motion of Stars in Declination')
plt.title('Plot of Proper Motion of Stars in Pleiades Cluster Data')
plt.show()

# Filtering out the stars that form a lump in proper motion space


# defining the sort function
def cutpm(dat, lowx, highx, lowy, highy, xaxis, yaxis):
    posx1 = np.where(dat[:, xaxis] < lowx)
    dat1 = np.delete(dat, posx1, 0)
    posx2 = np.where(dat1[:, xaxis] > highx)
    dat2 = np.delete(dat1, posx2, 0)
    posy1 = np.where(dat2[:, yaxis] < lowy)
    dat3 = np.delete(dat2, posy1, 0)
    posy2 = np.where(dat3[:, yaxis] < highy)
    dat4 = np.delete(dat3, posy2, 0)
    return dat4


# defining the cut data
cdatacut1 = cutpm(cdata, 17.5, 22, -48.5, -42, 7, 9)


# defining the function that gets the distances of stars from their parallax
def dist (p):
    d = 1 / (p / 1000)
    return d


# defining the distance data for the remaining stars
dists = dist(cdatacut1[:, 5])
#print(np.size(dists))
# plotting the histogram of stars vs distance from the cut data
plt.hist(dists, 50)
plt.xlabel('Distance to Stars [pc]')
plt.ylabel('Number of Stars')
plt.title('Histogram of Stars in Cluster with Proper Motions: \n RA: 17.5 - 22 mas/yr, Dec: -48.5 - -42 mas / yr')
plt.show()

# cutting down the data to make the histogram appear gaussian in shape and thus better


# defining the function to cut data based on distance/parallax
def cutpara(dat, low, high, axis, boundtype):
    if boundtype == 'p':
        pos1 = np.where(dat[:, axis] < low)
        dat1 = np.delete(dat, pos1, 0)
        pos2 = np.where(dat1[:, axis] > high)
        dat2 = np.delete(dat1, pos2, 0)
        return dat2
    elif boundtype == 'd':
        l = 1000 / high
        h = 1000 / low
        pos1 = np.where(dat[:, axis] < l)
        dat1 = np.delete(dat, pos1, 0)
        pos2 = np.where(dat1[:, axis] > h)
        dat2 = np.delete(dat1, pos2, 0)
        return dat2
    else:
        return 'Invalid'


# cutting the histogram data down
gauss_center = 250
sig2 = 0.49 * gauss_center
distlow = gauss_center - sig2
disthigh = gauss_center + sig2
cdatacut2 = cutpara(cdatacut1, distlow, disthigh, 5, 'd')

# plotting the revised histogram data
dists2 = dist(cdatacut2[:, 5])
print(np.size(dists2))
plt.hist(dists2, 30)
plt.xlabel('Distance to Stars [pc]')
plt.ylabel('Number of Stars')
plt.title('Histogram of Stars in Cluster after a 2 sigma cut')
plt.show()


# defining the magnitude function which calculates the absolute magnitudes of the luminosity of stars
def mag(dat, paraxis, magaxis):
    d = 1 / ((dat[:, paraxis]) / 1000)
    n = d / 10
    absmag = dat[:, magaxis] - 5 * np.log10(n)
    return absmag


# plotting the remaining stars on an HR diagram
colourcdc2 = cdatacut2[:, 14]
magcdc2 = mag(cdatacut2, 5, 12)

plt.plot(colourcdc2, magcdc2, 'o')
plt.ylim(20, -7)
plt.xlim(-0.5, 4)
plt.xlabel('B-R Colour')
plt.ylabel('Absolute Magnitude [M]')
plt.title('HR Diagram of Stars in Cluster')
plt.show()

# counting the number of white dwarves apparent on the HR diagram and returns the amount and the WD data
def countwd(dat, maglim, colourlim, magaxis, colouraxis):
    posmag = np.where(dat[:, magaxis] < maglim)
    dat1 = np.delete(dat, posmag, 0)
    poscol = np.where(dat1[:, colouraxis] < colourlim)
    n = np.size(poscol)
    return n


# getting the data for white dwarfs
def wddata(dat, maglim, colourlim, magaxis, colouraxis):
    posmag = np.where(dat[:, magaxis] < maglim)
    dat1 = np.delete(dat, posmag, 0)
    poscol = np.where(dat1[:, colouraxis] > colourlim)
    dat2 = np.delete(dat1, poscol, 0)
    return dat2


print('The numer of WDs in this HR diagram is ' + str(countwd(cdatacut2, 10, 0.7, 12, 14)))

# cutting excess stars and replotting the gaussian diagram
cdatacut3 = cutpara(cdatacut2, 50, 200, 5, 'd')


dists3 = dist(cdatacut3[:, 5])
print(np.size(dists3))
wds = wddata(cdatacut3, 10, 0.7, 12, 14)
distswd = dist(wds[:, 5])

plt.hist(dists3, 30)
plt.hist(distswd, 30)
plt.xlabel('Distance to Stars [pc]')
plt.ylabel('Number of Stars')
plt.title('Histogram of Stars in cluster after second cutting')
plt.show()


colourcdc3 = cdatacut3[:, 14]
magcdc3 = mag(cdatacut3, 5, 12)

magwds = mag(wds, 5, 12)
colourwds = wds[:, 14]

plt.plot(colourcdc3, magcdc3, 'o')
plt.ylim(17, -1)
plt.xlabel('B-R Colour')
plt.ylabel('Absolute Magnitude [M]')
plt.title('HR Diagram of Stars in Cluster')

plt.plot(colourwds, magwds, 'oy')
plt.show()

print('The numer of WDs in the revised set of data is ' + str(countwd(cdatacut3, 10, 0.7, 12, 14)))


# defining an averaging function to average extinctions while ignoring extinctions of 0 in value
def avg(dat, axis):
    data = dat[:, axis]
    count = 0
    val = 0
    for n in range(np.size(data)):
        if data[n] != 0:
            val = val + data[n]
            count += 1
        else:
            val = val + 0
            count = count + 0
    average = val / count
    return average


# Calculating the average extinction of the remaining stars in the cluster
avgext = avg(cdatacut3, 15)
print('the average extinction in g magnitude is '+ str(avgext))

avgbr = avg(cdatacut3, 16)
print('the average extinction in bp-rp magnitude is '+ str(avgbr))

# fixing and plotting the fixed colour data on an HR diagram
colourcdc3fix = cdatacut3[:, 14] - avgbr
magcdc3fix = mag(cdatacut3, 5, 12) - avgext

magwdsfix = mag(wds, 5, 12) - avgext
colourwdsfix = wds[:, 14] - avgbr

plt.plot(colourcdc3fix, magcdc3fix, 'o')
plt.ylim(20, -5)
plt.xlim(-0.5, 4)
plt.xlabel('B-R Colour')
plt.ylabel('Absolute Magnitude [M]')
plt.title('HR Diagram of Stars in Cluster After Extinction Correction')

plt.plot(colourwdsfix, magwdsfix, 'oy')
plt.show()


# plotting the best fitting isochrone atop the current HR diagram
"""
iso = 'is30.csv'
isochrone = np.loadtxt(iso, delimiter=',', comments='#', usecols=(23, 24, 25))
isomag = isochrone[:, 0]
isobprp = isochrone[:, 1] - isochrone[:, 2]

iso2 = 'is31.csv'
isochrone2 = np.loadtxt(iso2, delimiter=',', comments='#', usecols=(23, 24, 25))
isomag2= isochrone2[:, 0]
isobprp2 = isochrone2[:, 1] - isochrone2[:, 2]

iso3 = 'is32.csv'
isochrone3 = np.loadtxt(iso3, delimiter=',', comments='#', usecols=(23, 24, 25))
isomag3= isochrone3[:, 0]
isobprp3 = isochrone3[:, 1] - isochrone3[:, 2]

iso4 = 'is33.csv'
isochrone4 = np.loadtxt(iso4, delimiter=',', comments='#', usecols=(23, 24, 25))
isomag4= isochrone4[:, 0]
isobprp4 = isochrone4[:, 1] - isochrone4[:, 2]
"""

iso5 = 'is34.csv'
isochrone5 = np.loadtxt(iso5, delimiter=',', comments='#', usecols=(3, 23, 24, 25))
isomag5 = isochrone5[:, 1]
isobprp5 = isochrone5[:, 2] - isochrone5[:, 3]

plt.plot(colourcdc3fix, magcdc3fix, 'o')
# plt.plot(isobprp, isomag,  markersize=0.1)
# plt.plot(isobprp2, isomag2,  markersize=0.1)
# plt.plot(isobprp3, isomag3,  markersize=0.1)
# plt.plot(isobprp4, isomag4,  markersize=0.1)
plt.plot(isobprp5, isomag5,  markersize=0.1)
plt.ylim(25, -10)
plt.xlabel('B-R Colour')
plt.ylabel('Absolute Magnitude [M]')
plt.title('HR Diagram of Stars in Cluster including Isochrone')

plt.plot(colourwdsfix, magwdsfix, 'oy')
plt.show()

# defining the turnoff mass of the isochrone
massiso = isochrone5[:, 0]
turnoffmass = (massiso[96] + massiso[96] + massiso[96]) / 3
print('The turnoff mass is ' + str(turnoffmass) + 'solar masses')

# plotting the white dwarf cooling sequence of the cluster
wdcool = 'Table_Mass_07H.csv'
wdc = np.loadtxt(wdcool, delimiter=',', comments='#', usecols=(17, 18, 19))
wdcmag = wdc[:, 0]
wdcbprp = wdc[:, 1] - wdc[:, 2]

wdcool2 = 'Table_Mass_07He.csv'
wdc2 = np.loadtxt(wdcool2, delimiter=',', comments='#', usecols=(17, 18, 19))
wdcmag2 = wdc2[:, 0]
wdcbprp2 = wdc2[:, 1] - wdc2[:, 2]

wdcool3 = 'Table_Mass_12H.csv'
wdc3 = np.loadtxt(wdcool3, delimiter=',', comments='#', usecols=(17, 18, 19))
wdcmag3 = wdc3[:, 0]
wdcbprp3 = wdc3[:, 1] - wdc3[:, 2]

wdcool4 = 'Table_Mass_12He.csv'
wdc4 = np.loadtxt(wdcool4, delimiter=',', comments='#', usecols=(17, 18, 19))
wdcmag4 = wdc4[:, 0]
wdcbprp4 = wdc4[:, 1] - wdc4[:, 2]

plt.plot(colourcdc3fix, magcdc3fix, 'o')
plt.plot(wdcbprp, wdcmag,  markersize=0.1)
plt.plot(wdcbprp2, wdcmag2,  markersize=0.1)
plt.plot(wdcbprp3, wdcmag3,  markersize=0.1)
plt.plot(wdcbprp4, wdcmag4,  markersize=0.1)
plt.ylim(25, -10)
plt.xlabel('B-R Colour')
plt.ylabel('Absolute Magnitude [M]')
plt.title('HR Diagram of Stars in Cluster with WD Cooling Sequence')

plt.plot(colourwdsfix, magwdsfix, 'oy')
plt.show()

# Determining the main sequence lifetime of the white dwarfs
print(magwdsfix)
print(colourwdsfix)
age_cluster = 1.78E8  # years
#age_wd1 = 0  # years
#ms_lt1 = age_cluster - age_wd1
#age_wd2 = 0  # years
#ms_lt2 = age_cluster - age_wd2
age_wd3 = 98970000  # years
ms_lt3 = age_cluster - age_wd3
#print('The main sequence lifetime for the first white dwarf in this cluster is ' + str(ms_lt1) + ' years')
#print('The main sequence lifetime for the second white dwarf in this cluster is ' + str(ms_lt2) + ' years')
print('The main sequence lifetime for the third white dwarf in this cluster is ' + str(ms_lt3) + ' years')

# determining the mass of the white dwarf precursor
sun_mslt = 10E9  # years

#mass_prewd1 = (ms_lt1 / sun_mslt) ** (-2/5)
#print('The mass of the first WD precursor was ' + str(mass_prewd1) + ' solar masses')

#mass_prewd2 = (ms_lt2 / sun_mslt) ** (-2/5)
#print('The mass of the second WD precursor was ' + str(mass_prewd2) + ' solar masses')

mass_prewd3 = (ms_lt3 / sun_mslt) ** (-2/5)
print('The mass of the third WD precursor was ' + str(mass_prewd3) + ' solar masses')

# plotting the final vs initial masses of the white dwarfs in terms of solar masses
initial_mass = [mass_prewd3]
init_err = [0.4]
final_mass = [1.2]
final_err = [0.1]
plt.errorbar(initial_mass, final_mass, xerr=init_err, yerr=final_err, marker='s')
plt.xlabel('Initial Mass in solar masses [Msun}')
plt.ylabel('Final Mass in solar masses [Msun]')
plt.title('Initial-Final Mass Relation of White Dwarfs')
plt.show()

