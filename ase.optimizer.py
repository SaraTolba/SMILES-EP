from ase.calculators.dftb import Dftb
from ase.optimize import QuasiNewton,fire,gpmin,mdmin,LBFGS  
from ase.io import write,read
from ase.build import molecule
from rdkit import Chem
from rdkit.Chem import AllChem
from ase.optimize.minimahopping import MinimaHopping
from ase import *
from ase.optimize.basin import BasinHopping


def SmilesToMol(smiles):
    mol = Chem.MolFromSmiles(smiles)
    hmol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(hmol)
    AllChem.MMFFOptimizeMolecule(hmol)
    Chem.MolToMolFile(hmol,'trail.mol')
    return None
SmilesToMol(smiles='CNCc1cnccc1')
test = read('trail.mol') 
test.set_calculator(Dftb(label='h2o',
                         atoms=test,
                         Hamiltonian_MaxAngularMomentum_='',
                         Hamiltonian_MaxAngularMomentum_C='"p"',
                         Hamiltonian_MaxAngularMomentum_N='"p"',
                         Hamiltonian_MaxAngularMomentum_O='"p"',
                         Hamiltonian_MaxAngularMomentum_H='"s"',
                         ))
Quasi   =False
FIRE    =False
GPMin   =False
MDMin   =False
MinHop  =False
BasinHop=True
if Quasi: 
    dyn = QuasiNewton(test, trajectory='test.traj')
    dyn.run(fmax=0.001)
    write('test.final.xyz', test)
if FIRE:
    dyn = fire.FIRE(test, trajectory='test.traj')
    dyn.run()
    write('test.final.xyz', test)
if GPMin: 
    dyn = gpmin.gpmin.GPMin(test, trajectory='test.traj')
    dyn.run()
    write('test.final.xyz', test)
if MDMin:
    dyn = mdmin.MDMin(test, trajectory='test.traj')
    dyn.run()
    write('test.final.xyz', test)
if MinHop:
   opt = MinimaHopping(atoms=test)
   opt(totalsteps=20)
if BasinHop:
    bh = BasinHopping(atoms=test,         # the system to optimize
                      temperature=200, # 'temperature' to overcome barriers
                      dr=0.5,               # maximal stepwidth
                      optimizer=LBFGS,      # optimizer to find local minima
                      fmax=0.01,             # maximal force for the optimizer
                      )
