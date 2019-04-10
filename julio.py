from astroquery.gaia import Gaia
import astropy.coordinates as coord
import astropy.units as u
from astropy.coordinates import SkyCoord
from pylab import *

def get_time(V,s,A=1.1943):
	return (A*10**(V/5.)/s)**2

f = open('Downloads/temp_sample_gravity.txt','r')

teffs = []
bps   = []

lines = f.readlines()[1:]
ix = 0
tot = 0
ixf = 0
for line in lines:
	cos = line.split()
	ra1 = float(cos[1])
	dec1 = float(cos[2])
	ra2 = float(cos[12])
	dec2 = float(cos[13])
	bp1 = float(cos[-4])
	bp2 = float(cos[-2])

	#print ra1,dec1,ra2,dec2
	if dec1 < 30 and dec2 < 30:
		coord1 = SkyCoord(ra=ra1, dec=dec1, unit=(u.deg, u.deg), frame='icrs')
		coord2 = SkyCoord(ra=ra2, dec=dec2, unit=(u.deg, u.deg), frame='icrs')

		r1 = Gaia.query_object_async(coordinate=coord1, radius=u.Quantity(10./3600., u.deg))
		print len(r1)
		pmra = r1[0]['pmra']
		pmdec = r1[0]['pmdec']
		parallax = r1[0]['parallax']
		epmra = r1[0]['pmra_error']
		epmdec = r1[0]['pmdec_error']
		eparallax = r1[0]['parallax_error']	

		drop1,drop2 = False, False
		drop3,drop4 = False, False

		if len(r1)>1:
			for r in r1[1:]:
				if  (np.absolute(r['pmra']-pmra)/max(epmra,r['pmra_error']) < 5) and (np.absolute(r['pmdec']-pmdec)/max(epmdec,r['pmdec_error'])<5) and (np.absolute(r['parallax']-parallax)/eparallax < 5):
					drop1 = True
					break

		r2 = Gaia.query_object_async(coordinate=coord2, radius=u.Quantity(10./3600., u.deg))
		print len(r2)
		pmra = r2[0]['pmra']
		pmdec = r2[0]['pmdec']
		parallax = r2[0]['parallax']
		epmra = r2[0]['pmra_error']
		epmdec = r2[0]['pmdec_error']
		eparallax = r2[0]['parallax_error']	

		if len(r2)>1:
			for r in r2[1:]:
				if  (np.absolute(r['pmra']-pmra)/max(epmra,r['pmra_error']) < 5) and (np.absolute(r['pmdec']-pmdec)/max(epmdec,r['pmdec_error'])<5) and (np.absolute(r['parallax']-parallax)/eparallax < 5):
					drop2 = True
					break

		if r1['teff_val'][0]>7000:
			drop3 = True
		if r2['teff_val'][0]>7000:
			drop4 = True

		if drop1 or drop2 or drop3 or drop4:
			print 'OUT!'
		else:
			t1 = max(120.,get_time(bp1,20.))
			t2 = max(120.,get_time(bp2,20.))
			print ix
			print 'Teff1=',r1['teff_val'][0],'Teff2=',r2['teff_val'][0]
			print 'BP1=',bp1,'BP2=',bp2,'\n'
			print 'Texp1=',t1,'Texp2=',t2,'\n'
			tot += t1+t2


			plot([bp1,bp2],[r1['teff_val'][0],r2['teff_val'][0]],'k-')

			teffs.append(r1['teff_val'][0])
			teffs.append(r2['teff_val'][0])
			bps.append(bp1)
			bps.append(bp2)
			ixf += 1
			'\n'
		ix += 1

print 'final number of systems:', ixf
print 'total exposure time [hours]:', tot/3600.
xlabel('Vmag')
ylabel('Teff [K]')
scatter(bps,teffs,s=25,c='black',alpha=0.9)
show()