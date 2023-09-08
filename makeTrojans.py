
from spacerocks.simulation import Simulation
import pandas as pd
import argparse as ap
import rebound
import os

from astropy import units as u
import numpy as np

from itertools import groupby

from spacerocks import SpaceRock, Units
from astropy.time import Time

ref_epoch = Time('2023-08-11', format='iso', scale='tdb').jd

def has_adjacent_empty_cells(vec):
    return np.any(vec * np.roll(vec, 1))

def longest_contiguous_segment_1d(vec, target=True):
    return max(sum(g) for k, g in groupby(vec) if k == target)

def compute_longest_streak(mask, target=True):
    longest_duration = 0
    for row in mask.T:
        if has_adjacent_empty_cells(row):
            duration = longest_contiguous_segment_1d(row, target=target)
            if duration > longest_duration:
                longest_duration = duration
    return longest_duration

def number_of_segments(vec, target=True):
    return len([sum(g) for k, g in groupby(vec) if k == True])


def outdir(dirname):
    target = f'{dirname}'
    if os.path.isdir(target) is False:
        if os.path.exists(target):
            print(f'Error, bad output directory {target}')
        else:
            os.mkdir(target)
    print(f'Output dir = {target}')
    return target

def initial_population_clones(rock_to_clone, runid, N=10):
    name = rock_to_clone.name[0].replace(' ', '_')
    print(f'Making {N} clones of {name}...')
    clones = rock_to_clone.clone(N)
    output_dir = outdir(runid)
    clones.write_to(f'{output_dir}/initial_clones_{name}.rock')
    return clones

def initial_population_real(MPC_designations, runid):
    print(f'Retrieving {MPC_designations} from Horizons...')
    rocks = SpaceRock.from_horizons(MPC_designations)
    output_dir = outdir(runid)
    for name, rock in rocks.groupby('name'):
        pass   # The write currently crashes
        #rock.write_to(f'{output_dir}/initial_{name}.rock', compression=None)
    return rocks

def initial_population_synthetic(runid, outfile, N=5000, Ln=None):
#    ref_epoch = Time('2021-04-01', format='iso', scale='tdb').jd
    units = Units()
    units.timescale = 'tdb'
    units.timeformat = 'jd'
    units.mass = u.M_sun

    output_dir = outdir(runid)


    sim = Simulation(model='GIANTS', epoch=ref_epoch)
    _, planets, _ = sim.propagate(epochs=ref_epoch, exact_finish_time=1)
    jup = planets[planets.name == "Jupiter Barycenter"]
    mean_longitude_N = jup.mean_longitude.deg[0]

    rng = np.random.default_rng()

    arg = rng.uniform(0, 360, N)
    M = np.random.uniform(0, 360, N)
    # Choose node so that it lands in the libration region
    if Ln == "L4":
        node = mean_longitude_N - arg - M + rng.uniform(20, 110, size=N)
    elif Ln == "L5":
        node = mean_longitude_N - arg - M - rng.uniform(20, 110, size=N)
    else:  # make both if not specified
        sign = np.random.randint(0, 2, size=N) * 2 - 1
        node = mean_longitude_N - arg - M + sign * rng.uniform(20, 110, size=N)
    node = node % 360

    rocks = SpaceRock(a=rng.uniform(5.0, 5.4, N),
                      e=rng.uniform(0, 0.25, N),
                      inc=rng.uniform(0, 45, N),
                      arg=arg,
                      node=node,
                      M=M,
                      epoch=np.repeat(ref_epoch, N))
    rocks.write_to(f'{output_dir}/initial_rocks_{outfile}.rock')
    return rocks

def simulate(rocks, Nyears):
#    ref_epoch = Time('2021-04-01', format='iso', scale='tdb').jd
    model = 'GIANTS'
    units = Units()
    units.timescale = 'tdb'
    units.timeformat = 'jd'
    units.mass = u.M_sun
    simulation = Simulation(model=model, epoch=ref_epoch, units=units)

 #   print(rocks.name, rocks.epoch)
    simulation.add_spacerocks(rocks)
    simulation.integrator = 'whfast'
    simulation.dt = 40

    name_map = {rebound.hash(n).value: n for n in simulation.perturber_names + simulation.testparticle_names}

    def callback(sim):
        n_remaining = len(sim.remaining_testparticles)
#        print(f'Epoch = {(sim.t - sim.epoch[0].jd) / 365.25e06} Myr, remaining particles: {n_remaining}')
        for n in sim.remaining_testparticles:
            try:
                p = sim.particles[n]
                if (p.a < 4.9) or (p.a > 5.5):
                    #print(f'{name_map[p.hash.value]}, {(sim.t - sim.epoch[0].jd) / 365.25e06}')
                    sim.remove(hash=n)
                    sim.remaining_testparticles = np.delete(sim.remaining_testparticles,
                                                            np.where(sim.remaining_testparticles == n))
            except Exception as e:
                sim.remaining_testparticles = np.delete(sim.remaining_testparticles,
                                                        np.where(sim.remaining_testparticles == n))

    t0 = simulation.epoch.tdb.jd[0]
    tf = t0 + Nyears * 365.25
    epochs = np.arange(t0, tf, Nyears/1000 * 365.25)

    prop, planets, sim = simulation.propagate(epochs=epochs, units=units, callback=callback)

    return prop, planets, sim

def classify(prop, planets, sim):
    Nparticles = len(sim.remaining_testparticles)
    print(f'{Nparticles} remaining after integration')
    jup = planets[planets.name == "Jupiter Barycenter"]
    trojans = {}
    horseshoes = {}
    n = 0
    for name in sim.remaining_testparticles:
        if n % 10 == 0:
            print(f'Testing particle {n} of {Nparticles}...')

        clone = prop[prop.name == name]
        phi = (clone.mean_longitude.deg - jup.mean_longitude.deg) % 360

        t = (clone.epoch.jd - clone.epoch.jd.min()) / 365.25e6
        tbins = 30

        counts, _, _ = np.histogram2d(t, phi, range=[[t.min(), t.max()], [0, 360]], bins=(tbins, 30))
        mask = counts == 0
        streak = compute_longest_streak(mask)

        if streak == tbins:
            if phi.max() - phi.min() > 180:
                horseshoes[name] = clone
            else:
                trojans[name] = clone
        n += 1

    Ntrojans = len(trojans.items())
    Nhorseshoes = len(horseshoes.items())
    print(f"Found {Ntrojans} trojans and {Nhorseshoes} horseshoes.")

    return trojans, horseshoes

def saveobjects(runid, classified, outfile):
    output_dir = outdir(runid)
    a_trojans = []
    e_trojans = []
    inc_trojans = []
    arg_trojans = []
    node_trojans = []
    M_trojans = []
    epoch_trojans = []
    x_trojans = []
    y_trojans = []
    z_trojans = []
    vx_trojans = []
    vy_trojans = []
    vz_trojans = []
    trojans_name = []

    for name, rock in classified.items():
        rock.write_to(f'{output_dir}/{name}.rock', compression=None)
        trojans_a = rock.a[0].au
        a_trojans.append(trojans_a)

        trojans_e = rock.e[0]
        e_trojans.append(trojans_e)

        trojans_inc = rock.inc[0].deg
        inc_trojans.append(trojans_inc)

        trojans_arg = rock.arg[0].deg
        arg_trojans.append(trojans_arg)

        trojans_node = rock.node[0].deg
        node_trojans.append(trojans_node)

        trojans_M = rock.M[0].deg
        M_trojans.append(trojans_M)

        trojans_epoch = rock.epoch[0].jd
        epoch_trojans.append(trojans_epoch)

        trojans_name.append(name)

        trojans_x = rock.x[0].value
        trojans_y = rock.y[0].value
        trojans_z = rock.z[0].value
        trojans_vx = rock.vx[0].value
        trojans_vy = rock.vy[0].value
        trojans_vz = rock.vz[0].value
        x_trojans.append(trojans_x)
        y_trojans.append(trojans_y)
        z_trojans.append(trojans_z)
        vx_trojans.append(trojans_vx)
        vy_trojans.append(trojans_vy)
        vz_trojans.append(trojans_vz)

    units = Units()
    units.timescale = 'tdb'
    units.timeformat = 'jd'
    units.mass = u.M_sun

    sim = Simulation(model='GIANTS', epoch=epoch_trojans[0], units=units)
    epochs = epoch_trojans[0]
    _, planets, sim = sim.propagate(epochs=epochs, exact_finish_time=0)
    jup = planets[planets.name == "Jupiter Barycenter"]

    classified_dict = t= {'name':trojans_name, 'a':a_trojans, 'e':e_trojans, 'inc':inc_trojans, 'node':node_trojans,
                          'arg':arg_trojans, 'M':M_trojans, 'epoch':epoch_trojans,
                          'x':x_trojans, 'y':y_trojans, 'z':z_trojans, 'vx':vx_trojans, 'vy':vy_trojans, 'vz':vz_trojans}
    df_classified = pd.DataFrame.from_dict(classified_dict)

    Ln = []
    for ind, row in df_classified.iterrows():
        if row['y'] - jup.y.au > 0:
            Ln.append('L4')
        else:
            Ln.append('L5')
    df_classified['Ln'] = Ln

    df_classified.to_csv(f'{output_dir}/{outfile}.csv')
    return

def main():
    parser = ap.ArgumentParser(description='Configure Jupiter Trojan creation')
    parser.add_argument('-N', '--Nobjects',
                        help="Number of initial objects to integrate, eg. 5000", type=int, default=0)
    parser.add_argument('-y', '--years', required=True, help="Number of years, eg. 100,000,000", type=int)
    parser.add_argument('-o', '--output', required=True, help="Name of output file")
    parser.add_argument('-L', '--L', help="Type of Trojan to generate, = L4 or L5. If unspecified, makes both kinds.", type=str)
    parser.add_argument('-r', '--runid', help='Descriptor for this run', type=str, required=True)
    parser.add_argument('-d', '--designation', help='MPC designation of object to simulate', default=None, type=str)
    parser.add_argument('-c', '--clones', help='Number of clones to simulate (only if --designation is specified)', type=int)
    opts = parser.parse_args()
    N = opts.Nobjects
    Nyears = opts.years
    outfile = opts.output
    Ln = opts.L
    runid = opts.runid
    designation = opts.designation
    Nclones = opts.clones

    if N > 0:
        rocks = initial_population_synthetic(runid, outfile, N, Ln)
    else:
        if designation is not None and Nclones is None:
            rocks = initial_population_real(designation, runid)
        elif designation is not None and Nclones > 0:
            rock_to_clone = SpaceRock.from_horizons(designation)
            rocks = initial_population_clones(rock_to_clone, runid, Nclones)
        else:
            print('No rocks to simulate!')

    prop, planets, sim = simulate(rocks, Nyears)
    jup = planets[planets.name == "Jupiter Barycenter"]
    out_directory = outdir(runid)
    out = f'{out_directory}/Jupiter.rock'
    print(f'Writing {len(jup)} epochs for Jupiter to {out}')
    jup.write_to(out, compression=None)
    trojans, horseshoes = classify(prop, planets, sim)
    if len(trojans.items()) > 0:
        saveobjects(runid, trojans, f'trojans_{outfile}')
    if len(horseshoes.items()) == 0:
        return
    saveobjects(runid, horseshoes, f'horseshoes_{outfile}')



if __name__ == '__main__':
    main()

