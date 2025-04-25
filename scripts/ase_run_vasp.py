#! /usr/bin/env python
# adopted from: https://wiki.portal.chalmers.se/chemicalphysics/pmwiki.php/Theory/VASP-and-NEB-ASE

from ase.io import read, write
import os
from ase.neb import NEB
from ase.calculators.vasp import Vasp
from ase.optimize import FIRE



os.environ['VASP_PP_PATH'] = 'path/to/pseudopotentials'



def run(path):

    submitdir = os.getcwd() + f'/{path}/'
    print(submitdir)
    with open(f'{submitdir}/RUNNING', mode='a'): pass
    
    # you need to prepare traj_init.xyz migration trajectory
    # in submitdir
    images = read(f'{submitdir}/traj_init.xyz', index = ':')

    for i, image in enumerate(images):
        calc = Vasp(
                #restart = True,
                directory = submitdir + '{:02}'.format(i),
                encut=520,
                xc='PBE',
                gga='PE',
                kpts  = (1,1,1),
                gamma = True, # Gamma-centered (defaults to Monkhorst-Pack)
                ismear=0,
                nelm=100,
                nelmin=6,
                lreal='auto',
                ncore = 24,
                algo = 'fast',
                sigma = 0.1,
                ibrion=-1,
                ediffg=-0.1,  # forces
                ediff=1e-5,   # energy conv.
                prec='normal',
                lwave = True,
                nsw=0,        # don't use the VASP internal relaxation, only use ASE
                ispin=1,
                command = 'mpiexec -n 24 vasp_std',
                txt = 'out_neb_' + '{:02}'.format(i) + '.log',
                setups={'base': 'recommended', 'Li': ''}, # recommended PP, Li 2s1
                )
        image.calc = calc

    n_steps = 100
    print('Optimizing source')
    qn_source = FIRE(images[0],
            logfile=submitdir+'qn_source.log',
            trajectory=submitdir+'source.xyz',
            restart =submitdir+ 'qn_source.json',
            )
    qn_source.run(steps = n_steps, fmax = .1)
    
    print('Optimizing target')
    qn_target = FIRE(images[-1],
            logfile=submitdir+'qn_target.log',
            trajectory=submitdir+'target.xyz',
            restart = submitdir + 'qn_target.json',
            )
    qn_target.run(steps = n_steps, fmax = .1)


    neb = NEB(images,
            parallel=False,
            k = 5.0, 
            climb = True,
            method = 'improvedtangent'
            )

    qn = FIRE(neb,
            logfile=submitdir+'qn.log',
            trajectory=submitdir+'neb.xyz',
            restart = submitdir + 'qn.json'
    )

    print('Optimizing band')
    for i, _ in enumerate(qn.irun(fmax = 0.1, steps=n_steps)):
        forces = [abs(im.get_forces()).max() for im in qn.optimizable.neb.iterimages()]
        step = '{:03}'.format(i)
        print(f'step #{step}, |fmax|={max(forces)} eV/A')
        write(f'{submitdir}/band_optim_step_{i}.traj', qn.optimizable.neb.images)
    write(f'{submitdir}/traj_final.xyz', qn.optimizable.neb.images)

    n_images = len(images)
    for i in range(n_images):
        os.remove(f'{submitdir}/' + '{:02}'.format(i) + '/WAVECAR')
        os.remove(f'{submitdir}/' + '{:02}'.format(i) + '/CHGCAR')
        os.remove(f'{submitdir}/' + '{:02}'.format(i) + '/CHG')
    with open(f'{submitdir}/FINISHED', mode='a'): pass
    
    if abs(neb.get_forces()).max() < 0.1:
        with open(f'{submitdir}/CONVERGED', mode='a'): pass


# there are folders mp_id_source_target_offfset.neb
# each folder contains traj_init.xyz file
# e.g. mp-1211498_4_5_1_-1_-1.neb/traj_init.xyz
folders = [f for f in os.listdir() if '.neb' in f]
for folder in folders:
    try:
        if 'RUNNING' in os.listdir(folder):
            print('Calculation already exists in', folder)
            continue
        else:
            try:
                run(folder)
            except:
            #    raise
                print('failed')
                continue
    except KeyboardInterrupt:
        break
