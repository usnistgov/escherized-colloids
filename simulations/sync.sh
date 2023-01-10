#!/bin/bash
rsync -avz --exclude='Simulation.ipynb' --exclude='sync.sh' nam4@raritan.nist.gov:/home/nam4/waving_man/T=0.45/IH78/* ./
