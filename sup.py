# -*- coding: utf-8 -*-
"""
Created on Sun Oct  8 13:07:52 2023

@author: djele
"""

import os
import subprocess

# készítse el a képet is?
make_img = True

# szükséges programok elérési útja
vmd = "vmd" #"C://Program Files//VMD//vmd.exe"
povray = "povray"

# koordináták
# atomicnumbers 
# ide írd be ami kell, ha a "X atomic number not defined!" hibát kapod
symbols = { 1 : "H", 
            6 : "C",
            7 : "N",
            8 : "O",
           16 : "S",
           17 : "Cl"} 

min_int = 1.0 # minimum inteztiás a felhangokhoz és kombinációs sávokhoz

def read_mol(inp_fn):
    atoms = [] # vegyjelek
    std_coord = [] # koordiánaták
    with open(inp_fn,"r") as inp:
        lines = [ line.replace("\n","").lstrip().rstrip().lower() for line in inp.readlines()]
        
    # Anhramonikus vagy nem?
    anharm = any("anharmonic" in line for line in lines)

    # keressük meg az optimálás végét
    # for idx, line in enumerate(lines):
    #     if "stationary point found" in line: # ebből valamiért kettő van, egyes fájlokban
    #         idx_opt = idx
    #         break # álljunk meg az elsőnél, mert a második fura ???        
    # lines = lines[idx_opt+1:] # csak az ezután lévő részt használjuk

    # keressük standard orientation-t
    # feltevés - az utolsó a jó
    idxs = []
    idx = 0
    for line in lines:
        if "standard orientation:" in line:
            idxs.append(idx)
        idx += 1
    natoms = 0
    for line in lines[idxs[-1]+5:]:
        if  "----" in line:
            break
        split = line.split()
        try:
            atoms.append(symbols[int(split[1])])
        except ValueError:
            print(f"{split[1]} atomic number is not defined!")
        std_coord.append([float(split[3]),float(split[4]),float(split[5])])
        natoms += 1
    nfreqs = 3*natoms - 6
    ncombs= int((nfreqs*(nfreqs-1))/2)
    
    #mol = ase.Atoms("".join(atoms),positions=std_coord)
    #write(inp_fn.replace(".log",".png"),mol)
    
    # frekvenciák
    if anharm:
        modes = []
        freq_harms = []
        freq_anharms = []
        int_harms = []
        int_anharms = []
        idx_freqs = lines.index("fundamental bands")
        for idx in range(idx_freqs+3,idx_freqs+3+nfreqs):
            split = lines[idx].split()
            modes.append(split[0])
            freq_harms.append(float(split[1]))
            freq_anharms.append(float(split[2]))
            int_harms.append(float(split[3]))
            int_anharms.append(float(split[4]))
                
        over_modes = []
        over_freq_harms = []
        over_freq_anharms = []
        over_int_anharms = []
        over_idx_freqs = lines.index("overtones")
        for idx in range(over_idx_freqs+3,over_idx_freqs+3+nfreqs):
            split = lines[idx].split()
            if float(split[3]) >= min_int:
                over_modes.append(split[0])
                over_freq_harms.append(float(split[1]))
                over_freq_anharms.append(float(split[2]))
                over_int_anharms.append(float(split[3]))
        
        comb_modes = []
        comb_freq_harms = []
        comb_freq_anharms = []
        comb_int_anharms = []
        comb_idx_freqs = lines.index("combination bands")
        for idx in range(comb_idx_freqs+3,comb_idx_freqs+3+ncombs):
            split = lines[idx].split()
            if float(split[4]) >= min_int:
                comb_modes.append(split[0] + " " + split[1])
                comb_freq_harms.append(float(split[2]))
                comb_freq_anharms.append(float(split[3]))
                comb_int_anharms.append(float(split[4]))

    else:
        # frekvenciák
        freq_harms = []
        int_harms = []
        for line in lines:
            if "frequencies --" in line:
                if "low" in line:
                    continue
                split = line.replace("frequencies --","").split()
                for freq in split:
                    freq_harms.append(float(freq))
            elif "ir inten    --" in line:
                split = line.replace("ir inten    --","").split()
                for i in split:
                    int_harms.append(float(i))
        
        freq_harms = freq_harms[-1*nfreqs:]
        int_harms = int_harms[-1*nfreqs:]
                
    # keressük a "HF" és ZPVE energiát
    next_line = False
    e_hf_str = ""
    for line in lines:
        if next_line and e_hf_str != "":
            e_hf_str += str(line.split("\\")[0])
            E_HF = -1*abs(float(e_hf_str))
            next_line = False
        if "hf=" in line:
            for item in line.split("\\"):
                if "hf=" in item:
                    if not "\\" in line.split("hf=")[1]:
                        e_hf_str = str(item.split("=")[1])
                        next_line = True
                    else:
                        E_HF = float(item.split("=")[1])
                        
        if "zero-point correction=" in line:
            zpve = float(line.split("=")[1].split()[0])
    
    # ábrázoljuk a molekulánkat

    xyz_fn = "coord.xyz"
    pov_fn = "pov.pov"
    img_fn = inp_fn.replace(".log","")
    
    with open(xyz_fn, "w") as xyz_file:
        xyz_file.write(f"{natoms}\n")
        xyz_file.write("Molecule\n")
        for i in range(natoms):
            x, y, z = std_coord[i]
            element = atoms[i]
            xyz_file.write(f"{element} {x} {y} {z}\n")
    
    if make_img:
        print("-----------------------------------------------------------------")
        print("-----------------------------------------------------------------")
        print("Make The Image:")
        subprocess.run([vmd, "-e", "draw_mol.vmd"])
            
        subprocess.call([povray,f'-I{pov_fn}',f'-O{img_fn}',
                                "+W1240","+H1240","+D","+X","+A","+FN"])
        print("-----------------------------------------------------------------")
        print("-----------------------------------------------------------------")
    
    # írjuk ki az eredményeket Latex táblázatokba
    out_fn = "table.tex"
    with open(out_fn,"a") as out:
        # egy táblázat a koordinátáknak, képnek és energiáknak
        print(f"File name: {inp_fn}")
        title = input("Give the title: ")
        out.write("\\newpage \n")
        out.write("\\begin{center}\n")
        out.write("\\begin{longtable}{lrrrr}\n")
        line  = "\multicolumn{5}{c}{"
        line += title
        line += "} \\\\ \n"
        out.write(line)
        out.write("\\\\ \n")
        
        line = ""
        i = 0
        for elem, coord in zip(atoms,std_coord):
            line += "{0:} & ".format(elem)
            if i == 0:  
                line += "{0: 3.6f} & {1: 3.6f} & {2: 3.6f} & \n".format(coord[0],coord[1],coord[2])
                line += "ENERGY = {0: 5.6f}".format(E_HF)
                line += " Hartree \\\\ \n"
            elif i == 1:
                line += "{0: 3.6f} & {1: 3.6f} & {2: 3.6f} & \n".format(coord[0],coord[1],coord[2])
                line += "ZPVE = {0: 5.6f}".format(zpve)
                line += " Hartree \\\\ \n"
            elif i == 2:
                line += "{0: 3.6f} & {1: 3.6f} & {2: 3.6f} & \n".format(coord[0],coord[1],coord[2])
                line += "\multirow{"
                line += str(natoms)
                line += "}{*}{\includegraphics[width=7cm]{FiguresS/"
                line += inp_fn.replace(".log",".png")
                line += "}} \\\\ \n"
    
            else:
                line += "{0: 3.6f} & {1: 3.6f} & {2: 3.6f} & \\\\ \n".format(coord[0],coord[1],coord[2])
            i += 1
        out.write(line)
        
        out.write("\\end{longtable}\n")
        
        # frekvenciák kiírása egy "longtable"-be
        if anharm:
            out.write("\\begin{longtable}{crrrr} \n")
    
            out.write("\multicolumn{1}{c}{Mode(Quanta)} & ")
            out.write("\multicolumn{1}{c}{E(Harmonic)} & ")
            out.write("\multicolumn{1}{c}{I(Harmonic)} & ")
            out.write("\multicolumn{1}{c}{E(Anharmonic)} & ")
            out.write("\multicolumn{1}{c}{I(Anharmonic)} \\\\ \n")
    
            for mode, freq_harm, int_harm, freq_anharm, int_anharm in zip(modes,freq_harms,int_harms,freq_anharms,int_anharms):
                 line  = "{0:} & ".format(mode)
                 line += "{0: 5.3f} & ".format(freq_harm)
                 line += "{0: 5.3f} & ".format(int_harm)
                 line += "{0: 5.3f} & ".format(freq_anharm)
                 line += "{0: 5.3f} \\\\ \n".format(int_anharm)
                 out.write(line)
            
            out.write("\\\\ \n")
    
            for mode, freq_harm, freq_anharm, int_anharm in zip(over_modes,over_freq_harms,over_freq_anharms,over_int_anharms):
                 line  = "{0:} & ".format(mode)
                 line += "{0: 5.3f} & ".format(freq_harm)
                 line += " & "
                 line += "{0: 5.3f} & ".format(freq_anharm)
                 line += "{0: 5.3f} \\\\ \n".format(int_anharm)
                 out.write(line)
          
            out.write("\\\\ \n")
            for mode, freq_harm, freq_anharm, int_anharm in zip(comb_modes,comb_freq_harms,comb_freq_anharms,comb_int_anharms):
                 line  = "{0:} & ".format(mode)
                 line += "{0: 5.3f} & ".format(freq_harm)
                 line += " & "
                 line += "{0: 5.3f} & ".format(freq_anharm)
                 line += "{0: 5.3f} \\\\ \n".format(int_anharm)
                 out.write(line)
    
        else:
            out.write("\\begin{longtable}{rrr} \n")
    
            out.write("\multicolumn{1}{c}{Mode(Quanta)} & ")
            out.write("\multicolumn{1}{c}{E(Harmonic)} & ")
            out.write("\multicolumn{1}{c}{I(Harmonic)} \\\\ \n")
            i = 1
            for freq_harm, int_harm in zip(freq_harms,int_harms):
                line  = "{0: 3d}(1) & ".format(i)
                line += "{0: 5.3f} & ".format(freq_harm)
                line += "{0: 5.3f} \\\\ \n".format(int_harm)
                out.write(line)
                i += 1
    
        out.write("\\end{longtable}")
    
        out.write("\\end{center}\n")

# read_mol("./indene_anharm_Cs.log")

if os.path.isfile("table.tex"):
    os.remove("table.tex")

base = "."
for inp_fn in os.listdir(base):
    #inp_fn = os.path.join(base, fn)
    if os.path.isfile(inp_fn) and inp_fn.split(".")[-1] == "log":
        print(f"read info from: {inp_fn}")
        read_mol(inp_fn)
        
    
