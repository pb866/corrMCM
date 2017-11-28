"""
# Module C4CONO3CO

Correct MCMv3.3.1 mechanism for C4CONO3CO (and CO2N3CHO) species.

# Functions (public)
- both_rxn
"""
module C4CONO3CO

export both_rxn

using chk_mech


"""
    both_rxn(mech,beg_rxn,spc_end,newspc)

Correct species definitions and reactions for mechanism `mech` with both C4CONO3CO
and CO2N3CHO species and add degradation reactions for possible new products
as indicated in `newspc`.

Adjust parameters `spc_end` and `beg_rxn` indicating the line with the last
species definition and first reaction, respectively.

Returns mech, newspc, spc_end, and beg_rxn.
"""
function both_rxn(mech,beg_rxn,spc_end,newspc)

  # Loop over reactions and replace with corrections/delete obsolete reactions.
  for i = length(mech):-1:beg_rxn
    if contains(mech[i], "C4CONO3O2 = C4CONO3CO :")
       mech[i] = replace(mech[i],"C4CONO3O2 = C4CONO3CO :","C4CONO3O2 = CO2N3CHO :")
    elseif contains(mech[i], "C4CONO3CO + OH = C4CONO3O2") ||
      contains(mech[i], "C4CONO3CO + NO3 = C4CONO3O2")
      mech[i] = replace(mech[i],"C4CONO3O2","CO23C3CHO + NO2")
    elseif contains(mech[i], "C4CONO3CO =  C4CONO3O2 + CO + HO2 :")
       mech[i] = replace(mech[i],"C4CONO3CO =  C4CONO3O2 + CO + HO2 : 	J(34) 	;",
                 "C4CONO3CO = CH3CHO + HCOCO3 + NO2 :  J(34)+J(53) ;")
    elseif contains(mech[i], "C4CONO3CO = C4CO2O + NO2 : 	J(55) 	;")
      deleteat!(mech,i)
    end
  end

  # Add degradation reactions for possible new products
  # (only for species indicate in array of last parameter)
  mech,newspc,spc_end,beg_rxn = add_prod(mech,newspc,spc_end,beg_rxn,[1,2,5,6,7,8])

  # Return revised mechanism
  return mech,newspc,spc_end,beg_rxn
end #function COHM2_rxn

end #module C4CONO3CO
