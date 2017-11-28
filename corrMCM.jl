#!/usr/local/bin/julia


"""
# Module corrMCM

Correct duplicate species and their reactions in the MCMv3.3.1.
"""
module corrMCM

# Add path of self-made modules and load all modules
push!(LOAD_PATH,"/Applications/bin/data/jl.mod")
push!(LOAD_PATH,"~/Util/auxdata/jl.mod")
push!(LOAD_PATH,"./jl.mod")
using fhandle, chk_mech
import COHM2vCHOMOH, C4CONO3CO, NC3OO

# Add empty strings for missing arguments in ARGS
for i = 1:1-length(ARGS)  push!(ARGS,"")  end

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#

# Check existance of kpp file
kppfile = test_file(ARGS[1],default_dir="../..")

# Open temporary file for revised mechanism
tmpfile = splitext(kppfile)[1]*".tmp"
open(tmpfile,"w") do f

  ### Read kpp file
  mech = rdfil(kppfile)
  # Identify species in current mechanism for CHOMOH/COHM2 system
  nnew, newspc, spc_end, beg_rxn = index_mech(mech)

  #––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#

  ### Correction 1: C4CONO3CO & CO2N3CHO

  # Nothing to be done for CO2N3CHO, mechanism correct
  if nnew[1] > 0 && nnew[2] == 0
    println("Error from module C4CONO3COvCO2N3CHO:")
    println("No automatic corrections for mechanisms containing only C4CONO3CO currently available.")
    println("Please download mechanism from MCM website again adding CO2N3CHO to species list.")
    println("Script stopped"); exit()
  elseif nnew[1] > 0 && nnew[2] > 0
     mech,newspc,spc_end,beg_rxn = C4CONO3CO.both_rxn(mech,beg_rxn,spc_end,newspc)
  end


  #––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#

  ### Correction 2: The CHOMOHPAN/CHOMOHCO3/CHOMOHCO3H system

  # Replace mechanism according to case
  if nnew[3] > 0 && nnew[4] == 0
    mech,newspc,spc_end,beg_rxn = COHM2vCHOMOH.CHOMOH_rxn(mech,spc_end,beg_rxn,newspc)
  elseif nnew[4] > 0 && nnew[3] == 0
    # Mechanisms with only COHM2 species don't seem to exist
    println("Mechanism with only COHM2 species. Corrections not yet define.")
    println("Add routines to corrMCM.jl! Script stopped.")
    exit()
  elseif nnew[3] > 0 && nnew[4] > 0
    mech, spc_end, beg_rxn = COHM2vCHOMOH.del_spc(mech,beg_rxn,spc_end)
    mech,newspc,spc_end,beg_rxn = COHM2vCHOMOH.both_rxn(mech,beg_rxn,spc_end,newspc)
  end

  #––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#

  ### Correction 3: NOAOO & NC3OO

  # Nothing to be done for NOAOO and NOAOOA, mechanism correct
  if nnew[5] > 0 && nnew[6] > 0
    mech,spc_end,beg_rxn = NC3OO.del_spc(mech,spc_end,beg_rxn)
  end
  mech = NC3OO.replNC3OO(mech)

  #––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#

  ### Print revised mechanism to temp file
  for line in mech  println(f,line)  end
end #close tmpfile (f)
mv(tmpfile,kppfile,remove_destination=true)

end #module corrMCM
