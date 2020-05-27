import math
from datetime import datetime, timedelta
import time
import calendar


PI = math.pi
sin = math.sin
cos = math.cos
tan = math.tan
asin = math.asin
atan = math.atan2
acos = math.acos
rad = PI / 180.0

# date/time constants and conversions
dayMs = 1000 * 60 * 60 * 24
J1970 = 2440588
J2000 = 2451545
# calculations for sun times
J0 = 0.0009

# sun times configuration (angle, morning name, evening name)
times = [
    [-0.833, 'sunrise',       'sunset'      ],
    [  -0.3, 'sunriseEnd',    'sunsetStart' ],
    [    -6, 'dawn',          'dusk'        ],
    [   -12, 'nauticalDawn',  'nauticalDusk'],
    [   -18, 'nightEnd',      'night'       ],
    [     6, 'goldenHourEnd', 'goldenHour'  ]
]

e = rad * 23.4397  # obliquity of the Earth

toJulian = lambda dt: (calendar.timegm(dt.timetuple()) * 1000) / dayMs - 0.5 + J1970
fromJulian = lambda j: datetime.utcfromtimestamp(((j + 0.5 - J1970) * dayMs) / 1000.0)
toDays = lambda dt: toJulian(dt) - J2000

rightAscension = lambda l, b: atan(sin(l) * cos(e) - tan(b) * sin(e), cos(l))
declination = lambda l, b: asin(sin(b) * cos(e) + cos(b) * sin(e) * sin(l))
azimuth = lambda H, phi, dec: atan(sin(H), cos(H) * sin(phi) - tan(dec) * cos(phi))
altitude = lambda H, phi, dec: asin(sin(phi) * sin(dec) + cos(phi) * cos(dec) * cos(H))
siderealTime = lambda d, lw: rad * (280.16 + 360.9856235 * d) - lw

def astroRefraction(h):
    if h< 0:  # the following formula works for positive altitudes only.
        h = 0  # if h = -0.08901179 a div/0 would occur.
    # formula 16.4 of "Astronomical Algorithms" 2nd edition by Jean Meeus (Willmann-Bell, Richmond) 1998.
    # 1.02 / tan(h + 10.26 / (h + 5.10)) h in degrees, result in arc minutes -> converted to rad:
    return 0.0002967 / tan(h + 0.00312536 / (h + 0.08901179))

# general sun calculations


julianCycle = lambda d, lw: round(d - J0 - lw / (2 * PI))
approxTransit = lambda Ht, lw, n: J0 + (Ht + lw) / (2 * PI) + n
solarTransitJ = lambda ds, M, L: J2000 + ds + 0.0053 * sin(M) - 0.0069 * sin(2 * L)
hourAngle = lambda h, phi, d: acos((sin(h) - sin(phi) * sin(d)) / (cos(phi) * cos(d)))
observerAngle = lambda height: -2.076 * math.sqrt(height) / 60

solarMeanAnomaly = lambda d: rad * (357.5291 + 0.98560028 * d)


def hoursLater(date, h):
    return date + timedelta(hours=h)


def eclipticLongitude(M):
    C = rad * (1.9148 * sin(M) + 0.02 * sin(2 * M) + 0.0003 * sin(3 * M))  # equation of center
    P = rad * 102.9372  # perihelion of the Earth
    return M + C + P + PI


def sunCoords(d):
    M = solarMeanAnomaly(d)
    L = eclipticLongitude(M)
    return dict(dec=declination(L, 0), ra=rightAscension(L, 0))


# export functions
def getPosition(date, lat, lng):
    """calculates sun position for a given date and latitude/longitude"""
    lw = rad * -lng
    phi = rad * lat
    d = toDays(date)

    c = sunCoords(d)
    H = siderealTime(d, lw) - c["ra"]
    return dict(azimuth=azimuth(H, phi, c["dec"]), 
                altitude=altitude(H, phi, c["dec"]))


def addTime(angle, riseName, setName):
    times.append([angle, riseName, setName])


def getSetJ(h, lw, phi, dec, n, M, L):
    w = hourAngle(h, phi, dec)
    a = approxTransit(w, lw, n)
    return solarTransitJ(a, M, L)


# geocentric ecliptic coordinates of the moon
def moonCoords(d):
    L = rad * (218.316 + 13.176396 * d)
    M = rad * (134.963 + 13.064993 * d)
    F = rad * (93.272 + 13.229350 * d)

    l = L + rad * 6.289 * sin(M)
    b = rad * 5.128 * sin(F)
    dt = 385001 - 20905 * cos(M)

    return dict(ra=rightAscension(l, b), dec=declination(l, b), dist=dt)


def getMoonIllumination(date):
    d = toDays(date or datetime.now())
    s = sunCoords(d)
    m = moonCoords(d)

    # distance from Earth to Sun in km
    sdist = 149598000
    phi = acos(sin(s["dec"]) * sin(m["dec"]) + cos(s["dec"]) * cos(m["dec"]) * cos(s["ra"] - m["ra"]))
    inc = atan(sdist * sin(phi), m["dist"] - sdist * cos(phi))
    angle = atan(cos(s["dec"]) * sin(s["ra"] - m["ra"]),
                 sin(s["dec"]) * cos(m["dec"]) - cos(s["dec"]) * sin(m["dec"]) * cos(s["ra"] - m["ra"]))

    return dict(fraction=(1 + cos(inc)) / 2, phase=0.5 + 0.5 * inc * (-1 if angle < 0 else 1) / PI, angle=angle)


def getTimes(date, lat, lng, height):
    height = height or 0
    lw = rad * -lng
    phi = rad * lat

    dh = observerAngle(height)

    d = toDays(date)
    n = julianCycle(d, lw)
    ds = approxTransit(0, lw, n)

    M = solarMeanAnomaly(ds)
    L = eclipticLongitude(M)
    dec = declination(L, 0)

    Jnoon = solarTransitJ(ds, M, L)

    result = dict(solarNoon=fromJulian(Jnoon),
                nadir=fromJulian(Jnoon - 0.5)         
    )

    for i in range(0, len(times)):
        time = times[i]
        h0 = (time[0] + dh) * rad

        Jset = getSetJ(h0, lw, phi, dec, n, M, L)
        Jrise = Jnoon - (Jset - Jnoon)
        result[time[1]] = fromJulian(Jrise)
        result[time[2]] = fromJulian(Jset)

    return result


def getMoonPosition(date, lat, lng):
    lw = rad * -lng
    phi = rad * lat
    d = toDays(date)

    c = moonCoords(d)
    H = siderealTime(d, lw) - c["ra"]
    h = altitude(H, phi, c["dec"])
    pa = atan(sin(H), tan(phi) * cos(c['dec']) - sin(c['dec']) * cos(H))
    # altitude correction for refraction
    h = h + astroRefraction(h)

    return dict(azimuth=azimuth(H, phi, c["dec"]), altitude=h,
                 distance=c["dist"], parallacticAngle=pa)



def getMoonTimes(date, lat, lng):
    t = date.replace(hour=0, minute=0, second=0)

    hc = 0.133 * rad
    pos = getMoonPosition(t, lat, lng)
    h0 = pos["altitude"] - hc
    rise = 0
    sett = 0
    # go in 2-hour chunks, each time seeing if a 3-point quadratic curve crosses zero (which means rise or set)
    for i in range(1, 24, 2):
        h1 = getMoonPosition(hoursLater(t, i), lat, lng)["altitude"] - hc
        h2 = getMoonPosition(hoursLater(t, i + 1), lat, lng)["altitude"] - hc

        a = (h0 + h2) / 2 - h1
        b = (h2 - h0) / 2
        xe = -b / (2 * a)
        ye = (a * xe + b) * xe + h1
        d = b * b - 4 * a * h1
        roots = 0

        if d >= 0:
            dx = math.sqrt(d) / (abs(a) * 2)
            x1 = xe - dx
            x2 = xe + dx
            if abs(x1) <= 1:
                roots += 1
            if abs(x2) <= 1:
                roots += 1
            if x1 < -1:
                x1 = x2

        if roots == 1:
            if h0 < 0:
                rise = i + x1
            else:
                sett = i + x1

        elif roots == 2:
            rise = i + (x2 if ye < 0 else x1)
            sett = i + (x1 if ye < 0 else x2)

        if (rise and sett):
            break

        h0 = h2

    result = dict()

    if rise:
        result["rise"] = hoursLater(t, rise)
    if sett:
        result["set"] = hoursLater(t, sett)

    if (not rise and not sett):
        value = 'alwaysUp' if ye > 0 else 'alwaysDown'
        result[value] = True

    return result



if __name__ == '__main__':

    r = getTimes(datetime.strptime('2020-05-27', '%Y-%m-%d'), 44.4325, 26.1038889, 0)
    # r = getMoonTimes(datetime.strptime('2013-03-04', '%Y-%m-%d'), 50.5, 30.5, True)
    
    for k, v in r.items():
        print(k, v + timedelta(hours=3))
   