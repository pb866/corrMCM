"""
# Module COHM2vCHOMOH

Correct mechanism and substitute COHM2PAN/COHM2CO3H/COHM2CO3 species with
CHOMOHPAN/CHOMOHCO3H/CHOMOHCO3 species. Make adjustments to mechanism.

# Functions
## Public
- initCHOMOH
- initboth
- CHOMOH_rxn
- del_spc
- both_rxn

## Private
- subst_CHOMOH
"""
module COHM2vCHOMOH

export CHOMOH_rxn,
       del_spc,
       both_rxn

using chk_mech


"""
    CHOMOH_rxn(mech,spc_end,beg_rxn,newspc)

Correct species definitions and reactions for mechanism `mech` with only CHOMOH
species and add degradation reactions for new products as indicated in `newspc`.

Adjust parameters `spc_end` and `beg_rxn` indicating the line with the last
species definition and first reaction, respectively.

Returns mech, newspc, spc_end, and beg_rxn.
"""
function CHOMOH_rxn(mech,spc_end,beg_rxn,newspc)

  # substitute mechanism with corrections
  for i = beg_rxn:length(mech)
    mech[i] = subst_CHOMOH(mech[i])
  end
  # Define species, which to look for
  chk = Int64[]; for i = 3:8  push!(chk,i)  end
  # Add degradation reactions for possible new products
  mech,newspc,spc_end,beg_rxn = add_prod(mech,newspc,spc_end,beg_rxn,chk)

  # Return revised mechanism and indices
  return mech,newspc,spc_end,beg_rxn
end #function CHOMOH_rxn


"""
    del_spc(mech,beg_rxn,spc_end)

Delete species definitions and concentrations in the RO2 summation of COHM2 species
(except COHM2CO2H) in `mech`, which hold both COHM2 and CHOMOH species starting
at the first reaction with index `beg_rxn`. Adjust parameter `spc_end` with index
of last species definition.

Returns mech, spc_end, and beg_rxn.
"""
function del_spc(mech,beg_rxn,spc_end)

  # Loop over species definitions and RO2 summation
  # (in reverse order to not jump over lines, when deleting species definitions)
  for i = beg_rxn:-1:1
    # Delete species definitions of COHM2 species
    if ismatch(r"COHM2(?!CO2H).*IGNORE",mech[i])
      deleteat!(mech,i); beg_rxn -= 1; spc_end -= 1
    end
    # Delete COHM2 species in RO2 summation (consider several possible print forms)
    if ismatch(r" \+ C\(ind_COHM2[A-Za-z0-9 ]*\)",mech[i])
      mech[i] = replace(mech[i],r" \+ C\(ind_COHM2[A-Za-z0-9 ]*\)","")
      beg_rxn -= 1
    elseif ismatch(r"C\(ind_COHM2[A-Za-z0-9 ]*\) \+ ",mech[i])
      mech[i] = replace(mech[i],r"C\(ind_COHM2[A-Za-z0-9 ]*\) \+ ","")
      beg_rxn -= 1
    elseif ismatch(r"C\(ind_COHM2[A-Za-z0-9 ]*\)",mech[i])
      mech[i] = replace(mech[i],r"C\(ind_COHM2[A-Za-z0-9 ]*\) \+ ","")
      print("Warning! COHM2 RO2 detected on line $i ")
      println("in RO2 summation, but not automatically removed.")
    end
  end

  # Return revised mechanism and index for start of reactions
  return mech, spc_end, beg_rxn
end #function del_spc


"""
    both_rxn(mech,beg_rxn,spc_end,newspc)

Correct species definitions and reactions for mechanism `mech` with both CHOMOH
and COHM2 species and add degradation reactions for possible new products
as indicated in `newspc`.

Adjust parameters `spc_end` and `beg_rxn` indicating the line with the last
species definition and first reaction, respectively.

Returns mech, newspc, spc_end, and beg_rxn.
"""
function both_rxn(mech,beg_rxn,spc_end,newspc)

  # Loop over reactions
  for i = length(mech):-1:beg_rxn
    # Replace COHM2 with CHOMOH species
    if contains(mech[i], "MACRNBCO3H + OH = COHM2CO3H + NO2") ||
       contains(mech[i], "MACRNBPAN + OH = COHM2PAN + NO2") ||
       contains(mech[i], "COHM2PAN = COHM2CO3 + NO2 :") ||
       contains(mech[i], "COHM2CO3H + OH = COHM2CO3 :") ||
       contains(mech[i], "COHM2CO3 + HO2 = COHM2CO3H :") ||
       contains(mech[i], "COHM2CO3 + HO2 = COHM2CO2H + O3 :") ||
       contains(mech[i], "COHM2CO3 = COHM2CO2H :")
       mech[i] = replace(mech[i],r"COHM2(?!CO2H)","CHOMOH")

    # Delete duplicate COHM2 species
    elseif contains(mech[i], "COHM2CO3 + NO2 = COHM2PAN") ||
       contains(mech[i], "CHOMOHPAN = CHOMOHCO3 + NO2 : 	KBPAN 	;") ||
       contains(mech[i], "COHM2PAN + OH = GLYOX + NO3 :") ||
       contains(mech[i], "CHOMOHCO3H + OH = CHOMOHCO3 : 	6.99D-11 	;") ||
       contains(mech[i], "COHM2CO3H = GLYOX + HO2 + OH :") ||
       contains(mech[i], "COHM2CO3H = HCOCO3H + CO + HO2 :") ||
       contains(mech[i], "COHM2CO3 + HO2 = GLYOX + HO2 + OH :") ||
       contains(mech[i], "CHOMOHCO3 + HO2 = CHOMOHCO3H : 	KAPHO2*0.56 	;") ||
       contains(mech[i], "COHM2CO3 + NO = GLYOX + HO2 + NO2 :") ||
       contains(mech[i], "COHM2CO3 + NO3 = GLYOX + HO2 + NO2 :") ||
       contains(mech[i], "COHM2CO3 = GLYOX + HO2 :")
       deleteat!(mech,i)

    # Add branching ratio to RCO3 main channel
    elseif contains(mech[i], "CHOMOHCO3 = MGLYOX + HO2 :")
       mech[i] = replace(mech[i],r"(?<!0\.7\*)RO2","0.7*RO2")
    end
  end

  # Define species, which to look for
  chk = Int64[]; for i = 3:8  push!(chk,i)  end
  # Add degradation reactions for possible new products
  mech,newspc,spc_end,beg_rxn = add_prod(mech,newspc,spc_end,beg_rxn,chk)

  # Return revised mechanism
  return mech,newspc,spc_end,beg_rxn
end #function COHM2_rxn


"""
    subst_CHOMOH(line)

Substitute `line` with corrections of the revised mechanism.
"""
function subst_CHOMOH(line)

  if contains(line,"CHOMOHPAN = CHOMOHCO3 + NO2")
    line = replace(line,r"KBPAN(?!\+J\(17\))","KBPAN+J(17)")
  elseif contains(line,"CHOMOHCO3H + OH = CHOMOHCO3 :")
    line = replace(line,"6.99D-11","2.47D-11")
  elseif contains(line,"CHOMOHCO3 + HO2 = CHOMOHCO3H : 	KAPHO2*0.56")
    rxn = match(r"\{(.*?)\.",line)[1]
    line = replace(line,".}","a}",1)
    if !contains(line,"0.56")
      line = replace(line,"0.56","0.41")*"\n{"*rxn*
             "b} 	 CHOMOHCO3 + HO2 = COHM2CO2H + O3 : 	KAPHO2*0.15 	;"
    end
  elseif contains(line,"CHOMOHCO3 = MGLYOX + HO2")
    rxn = match(r"\{(.*?)\.",line)[1]
    if !contains(line,"0.7*RO2")
      line = replace(line,".}","a}",1)
      line = replace(line,"RO2","0.7*RO2")*"\n{"*rxn*
             "b} 	 CHOMOHCO3 = COHM2CO2H : 	1.00D-11*0.3*RO2 	;"
    end
  end

  return line
end #function subst_CHOMOH

end #module COHM2vCHOMOH
