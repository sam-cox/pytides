"""
This module is part of the Python package Pytides.

The MIT License (MIT)

Copyright (c) 2013 Sam Cox where applicable.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

from collections import OrderedDict, Iterable
from itertools import izip, takewhile, count, ifilter
from datetime import datetime, timedelta
import numpy as np
from scipy.optimize import leastsq, fsolve
from astro import astro
import constituent

d2r, r2d = np.pi/180.0, 180.0/np.pi

class Tide(object):
    dtype = np.dtype([
        ('constituent', object),
        ('amplitude', float),
        ('phase', float)])

    def __init__(self, constituents = None, amplitudes = None, phases = None,
                 model = None, radians = False):
        """Initialise a tidal model. """
        if None not in [constituents, amplitudes, phases]:
            if len(constituents) == len(amplitudes) == len(phases):
                model = np.zeros(len(phases), dtype=Tide.dtype)
                model['constituent'] = np.array(constituents)
                model['amplitude'] = np.array(amplitudes)
                model['phase'] = np.array(phases)
            else:
                raise ValueError("Constituents, amplitudes and phases should all be arrays of equal length.")
        elif model is not None:
            if not model.dtype == Tide.dtype:
                raise ValueError("Model must be a numpy array with dtype == Tide.dtype")
        else:
            raise ValueError("Must be initialised with constituents, amplitudes and phases; or a model.")
        if radians:
            model['phase'] = r2d*model['phase']
        self.model = model[:]
        self.normalize()
    
    def prepare(self, *args, **kwargs):
        return Tide._prepare(self.model['constituent'], *args, **kwargs)

    @staticmethod 
    def _prepare(constituents, t0, t = None, radians = True):
        #The equilibrium argument is constant and taken at the beginning of the time series (t0).
        #The speed of the equilibrium argument changes very slowly, so again we take it to be constant
        #over any length of data. The node factors change more rapidly.
        if isinstance(t0, Iterable):
            t0 = t0[0]
        if t is None:
            t = [t0]
        if not isinstance(t, Iterable):
            t = [t]
        a0 = astro(t0)
        a = [astro(t_i) for t_i in t]
        
        #For convenience give u, V0 (but not speed!) in [0, 360)
        V0 = np.array([c.V(a0) for c in constituents])[:, np.newaxis]
        speed = np.array([c.speed(a0) for c in constituents])[:, np.newaxis]
        u = [np.mod(np.array([c.u(a_i) for c in constituents])[:, np.newaxis], 360.0)
             for a_i in a]
        f = [np.mod(np.array([c.f(a_i) for c in constituents])[:, np.newaxis], 360.0)
             for a_i in a]
        
        if radians:
            speed = d2r*speed
            V0 = d2r*V0
            u = [d2r*each for each in u]
        return speed, u, f, V0
         
    def at(self, t):
        if isinstance(t, Iterable):
            #We'll do this a bit more sensibly soon
            return np.ndarray([self.at(ti) for ti in t])
        speed,[u],[f],V0 = self.prepare(t, radians=True)
        t = np.zeros(1)#self.hours(t)
        amplitudes = self.model['amplitude'][:, np.newaxis]
        phases = d2r*self.model['phase'][:, np.newaxis]
        return Tide._tidal_series(t, amplitudes, phases, speed,u,f,V0)
    
    def highs(self, args):
        for t in ifilter(lambda e: e[2] == 'H', self.extrema(*args)):
            yield t
            
    def lows(self, args):
        for t in ifilter(lambda e: e[2] == 'H', self.extrema(*args)):
            yield t
    
    def form_number(self):
        k1, o1, m2, s2 = (np.extract(self.model['constituent'] == c, self.model['amplitude']) for c
                          in [constituent._K1, constituent._O1, constituent._M2, constituent._S2])
        return (k1+o1)/(m2+s2)
    
    def classify(self):
        form = self.form_number()
        if 0 <= form <= 0.25:
            return 'semidiurnal'
        elif 0.25 < form <= 1.5:
            return 'mixed (semidiurnal)'
        elif 1.5 < form <= 3.0:
            return 'mixed (diurnal)'
        else:
            return 'diurnal'
        
    def extrema(self, t0, t1 = None, partition = 2400.0):
        if t1:
            #yield from in python 3.4
            for e in takewhile(lambda t: t[0] < t1, self.extrema(t0)):
                yield e
        else:
            delta = np.amin([90.0 / c.speed(astro(t0)) for c in self.model['constituent'] if not c.speed(astro(t0)) == 0])
            #We search for stationary points from offset hours before t0 to ensure we find any which
            #might occur very soon after t0.
            offset = 24.0
            amplitude = self.model['amplitude'][:, np.newaxis]
            phase = d2r*self.model['phase'][:, np.newaxis]
            partitions = (Tide._times(t0, i*partition) for i in count()), (Tide._times(t0, i*partition) for i in count(1))
            #We'll overestimate to be on the safe side; values outside (start,end) won't get yielded.
            interval_count = int(np.ceil((partition + offset) / delta)) + 1
            for start, end in izip(*partitions):
                speed, [u], [f], V0 = self.prepare(start, Tide._times(start, 0.5*partition))
                # These derivatives don't include the time dependence of u or f, but these change slowly.
                def d(t):
                    return np.sum(-speed*amplitude*f*np.sin(speed*t + (V0 + u) - phase), axis=0)
                def d2(t):
                    return np.sum(-speed**2.0 * amplitude*f*np.cos(speed*t + (V0 + u) - phase), axis=0)
                #We'll overestimate to be on the safe side; values outside (start,end) won't get yielded.
                intervals = (delta * i -offset for i in range(interval_count)), (delta*(i+1) - offset for i in range(interval_count))
                for a, b in izip(*intervals):
                    if d(a)*d(b) < 0:
                        extrema = fsolve(d, (a + b) / 2.0, fprime = d2)[0]
                        time = Tide._times(start, extrema)
                        height = self.at(time)
                        type = 'H' if d2(extrema) < 0 else 'L'
                        if start < time < end:
                            yield (time, height, type)
    
    @staticmethod
    def _hours(t0, t):
        if not isinstance(t, Iterable):
            return Tide._hours(t0, [t])[0]
        elif isinstance(t[0], datetime):
            return np.array([(ti-t0).total_seconds() / 3600.0 for ti in t])
        else:
            return t
        
    @staticmethod
    def _partition(hours, partition = 3600.0):
        partition = float(partition)
        relative = hours - hours[0]
        total_partitions = np.ceil(relative[-1] / partition + 10*np.finfo(np.float).eps).astype('int')
        return [hours[np.floor(np.divide(relative, partition)) == i] for i in range(total_partitions)]        
    
    @staticmethod
    def _times(t0, hours):
        if not isinstance(hours, Iterable):
            return Tide._times(t0, [hours])[0]
        elif not isinstance(hours[0], datetime):
            return np.array([t0 + timedelta(hours=h) for h in hours])
        else:
            return np.array(hours)
    ##Return the modelled heights for given times
    ##t is of shape (nt,); all other arguments are of shape (nc, 1)
    @staticmethod
    def _tidal_series(t, amplitude, phase, speed, u, f, V0):     
        return np.sum(amplitude*f*np.cos(speed*t + (V0 + u) - phase), axis=0)
    
    #By convention amplitudes, phases are positive
    def normalize(self):
        for i, (_, amplitude, phase) in enumerate(self.model):
            if amplitude < 0:
                self.model['amplitude'][i] = -amplitude
                self.model['phase'][i] = phase + 180.0    
            self.model['phase'][i] = np.mod(self.model['phase'][i], 360.0) 

    @classmethod
    def decompose(cls, heights, t = None, t0 = None, interval=None, constituents=constituent.noaa, initial = None, n_period = 2, callback = None, full_output = False):
        if t is not None:
            if isinstance(t[0], datetime):
                hours = Tide._hours(t[0], t)
                t0 = t[0]
            elif t0 is not None:
                hours = t
            else:
                raise ValueError("t can be an array of datetimes, or an array of hours since t0 in which case t0 must be specified.")       
        elif None not in [t0, interval]:
            hours = np.arange(len(heights)) * interval
        else:
            raise ValueError("Must provide t(datetimes), or t(hours) and t0(datetime), or interval(hours) and t0(datetime) so that each height can be identified with an instant in time.")

        #Remove duplicate constituents (those which travel at exactly the same speed, irrespective of phase)
        constituents = list(OrderedDict.fromkeys(constituents))

        #No need for least squares to find the mean water level constituent z0, work relative to mean
        constituents = [c for c in constituents if not c == constituent._Zo]
        z0 = np.mean(heights)
        heights = heights - z0

        #Only analyse frequencies which complete at least n_period cycles over the data period.
        constituents = [c for c in constituents if 360.0 * n_period < hours[-1] * c.speed(astro(t0))]
        n = len(constituents)
        
        sort = np.argsort(hours)
        hours = hours[sort]
        heights = heights[sort]
        
        #We partition our time/height data into intervals over which we consider the values of
        #u and f to assume a constant value (that is, their true value at the midpoint of the interval).
        #Constituent speeds change much more slowly than the node factors, so we will consider these
        #constant and equal to their speed at t0, regardless of the length of the time series.
        
        partition = 240.0
        t = Tide._partition(hours, partition)
        times = Tide._times(t0, [(i + 0.5)*partition for i in range(len(t))])
        speed, u, f, V0 = Tide._prepare(constituents, t0, times, radians = True)
        
        #Residual to be minimised by variation of parameters (amplitudes, phases)
        def residual(hp):
            H, p = hp[:n, np.newaxis], hp[n:, np.newaxis]    
            s = np.concatenate([Tide._tidal_series(t_i, H, p, speed, u_i, f_i, V0) for t_i, u_i, f_i in izip(t, u, f)])
            res = heights - s
            if callback:
                callback(res)
            return res

        ##Analytic Jacobian of the residual - this makes solving significantly faster than just using
        ##gradient approximation, especially with many measurements / constituents.
        def D_residual(hp):
            H, p = hp[:n, np.newaxis], hp[n:, np.newaxis]
            ds_dH = np.concatenate([f_i*np.cos(speed*t_i+u_i+V0-p) for t_i, u_i, f_i in izip(t, u, f)], axis = 1)
            ds_dp = np.concatenate([H*f_i*np.sin(speed*t_i+u_i+V0-p) for t_i, u_i, f_i in izip(t, u, f)], axis = 1)
            return np.append(-ds_dH, -ds_dp, axis=0)

        #Initial guess for solver, haven't done any analysis on this since the solver seems to converge well regardless of the initial guess
        #We do however scale the initial amplitude guess with some measure of the variation
        amplitudes = np.ones(n) * (np.sqrt(np.dot(heights, heights)) / len(heights))
        phases = np.ones(n)
        
        if initial:
            for (c0, amplitude, phase) in initial.model:
                for i, c in enumerate(constituents):
                    if c0 == c:
                        amplitudes[i] = amplitude
                        phases[i] = d2r*phase

        initial = np.append(amplitudes, phases)

        lsq = leastsq(residual, initial, Dfun=D_residual, col_deriv=True, ftol=1e-7) # 

        model = np.zeros(1+n, dtype=cls.dtype)
        model[0] = (constituent._Zo, z0, 0)
        model[1:]['constituent'] = constituents[:]
        model[1:]['amplitude'] = lsq[0][:n]
        model[1:]['phase'] = lsq[0][n:]

        if full_output:
            return cls(model = model, radians = True), lsq
        return cls(model = model, radians = True)
