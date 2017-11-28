"""
# Module COHM2vCHOMOH

Find species to correct in current mechanism and missing degradation reactions
for new products.

# Functions (public)
- index_mech
- add_prod
"""
module chk_mech

export index_mech,
       add_prod


"""
index_mech(mech)

Find number of occurances `nnew` of species with corrections in `mech` and
mark those species and their products, if present in the current mechanism,
in `newspc`. Find indices for end of species definitions (`spc_end`) and
beginning of reactions (`beg_rxn`) in `mech`.

Returns `nnew`, `newspc`, `spc_end`, and `beg_rxn`.
"""
function index_mech(mech)

  # Initialise
  newspc = falses(10)
  nnew = zeros(6)

  # Loop over kpp file
  for line in mech
    # Count occurances of species, where corrections are necessary
    nnew[1] += length(matchall(r"C4CONO3CO",line))
    nnew[2] += length(matchall(r"CO2N3CHO",line))
    nnew[3] += length(matchall(r"CHOMOH",line))
    nnew[4] += length(matchall(r"COHM2(?!(CO2H))",line))
    nnew[5] += length(matchall(r"NOAOO",line))
    nnew[6] += length(matchall(r"NC3OO",line))
    # Find products in kpp file from new reactions
    # Initialise flags for each new species
    # newspc[COHM2CO2H,GLYOX,HCOCO,HCOCO3,HCOCO2H,HCOCO3H]
    if contains(line," CO23C3CHO ") newspc[1] = true
    elseif contains(line," CH3CHO ") newspc[2] = true
    elseif contains(line," COHM2CO2H ") newspc[3] = true
    elseif contains(line," GLYOX ") newspc[4] = true
    elseif contains(line," HCOCO ") newspc[5] = true
    elseif contains(line," HCOCO3 ") newspc[6] = true
    elseif contains(line," HCOCO2H ") newspc[7] = true
    elseif contains(line," HCOCO3H ") newspc[8] = true
    elseif contains(line," CH3CO3") newspc[9] = true
    elseif contains(line," HCOCH2O2") newspc[10] = true
    end
  end

  # Find end of species definitions and beginning of reactions
  spc_end = find([contains(line,"{ Peroxy radicals. }") for line in mech])[1]
  beg_rxn = find([contains(line,"#EQUATIONS") for line in mech])[1] + 1

  return nnew, newspc, spc_end, beg_rxn
end #function sel_CHOMOH_case



"""
    add_prod(mech,newspc,spc_end,beg_rxn,chk::Array{Int64})

Add degradation reactions of possible new products to revised mechanism `mech`.

New products are flagged in `newspc`. The routine only checks for those species
in `newspc` whose index is given in `chk`. Furthermore the parameters for the
end of the species definitions and the beginning of reactions, `spc_end` and
`beg_rxn`, are adjusted, if new species are inserted.

Returns mech, newspc, spc_end, and beg_rxn.
"""
function add_prod(mech,newspc,spc_end,beg_rxn,chk::Array{Int64})

  # Find last reaction number
  for i = length(mech):-1:1
    try mech[i][1] == '{'
      rxn = parse(Int64,match(r"\{(.*?)\.",mech[end])[1])
      break
    end
  end

  #––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#

  # Products from correction 1:

  # Add reactions for CO23C3CHO, if currently not in the mechanism
  if !newspc[1] && any(x->x==1,chk)
    rxn += 1
    push!(mech, "{$rxn.} 	 CO23C3CHO = CH3CO3 + CO + CO + HO2 : 	J(34) 	;")
    rxn += 1
    push!(mech, "{$rxn.} 	 CO23C3CHO = CH3CO3 + HCOCO : 	J(35) 	;")
    rxn += 1
    push!(mech, "{$rxn.} 	 NO3 + CO23C3CHO = CH3CO3 + CO + CO + HNO3 : 	KNO3AL*4.0 	;")
    rxn += 1
    push!(mech, "{$rxn.} 	 OH + CO23C3CHO = CH3CO3 + CO + CO : 	1.23D-11 	;")
    # Set species and new products true
    newspc[1] = true; newspc[5:8] = true
    # Add species definitions
    insert!(mech,spc_end,"CO23C3CHO = IGNORE ;"); spc_end += 1; beg_rxn += 1
    # Warn about missing products, currently not checked:
    if !newspc[9]
      println("Warning! No degradation for 'CH3CO3' described in current mechanism.")
      print("Download new mechanism from MCM website or adjust ")
      println("function add_prod in module chk_mech.")
    end
  end
  # Add reactions for CH3CHO, if currently not in the mechanism
  if !newspc[2] && any(x->x==2,chk)
    rxn += 1
    push!(mech, "{$rxn.} 	 CH3CHO = CH3O2 + HO2 + CO : 	J(13) 	;")
    rxn += 1
    push!(mech, "{$rxn.} 	 NO3 + CH3CHO = HNO3 + CH3CO3 : 	1.4D-12*EXP(-1860/TEMP) 	;")
    rxn += 1
    push!(mech, "{$rxn.} 	 OH + CH3CHO = CH3CO3 : 	4.7D-12*EXP(345/TEMP)*0.95 	;")
    rxn += 1
    push!(mech, "{$rxn.} 	 OH + CH3CHO = HCOCH2O2 : 	4.7D-12*EXP(345/TEMP)*0.05 	;")
    rxn += 1
    push!(mech, "{$rxn.} 	 ")
    # Set species and new products true
    newspc[[3,4,7]] = true
    # Add species definitions
    insert!(mech,spc_end,"CH3CHO = IGNORE ;"); spc_end += 1; beg_rxn += 1
    # Warn about missing products, currently not checked:
    if !newspc[9]
      println("Warning! No degradation for 'CH3CO3' described in current mechanism.")
      print("Download new mechanism from MCM website or adjust ")
      println("function add_prod in module chk_mech.")
    end
    if !newspc[10]
      println("Warning! No degradation for 'HCOCH2O2' described in current mechanism.")
      print("Download new mechanism from MCM website or adjust ")
      println("function add_prod in module chk_mech.")
    end
  end

  #––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#

  # Products from correction 2:

  # Add reactions for COHM2CO2H, if currently not in the mechanism
  if !newspc[3] && any(x->x==3,chk)
    rxn += 1
    push!(mech, "{$rxn.} 	 COHM2CO2H + OH = GLYOX + HO2 : 	2.16D-11 	;")
    rxn += 1
    push!(mech, "{$rxn.} 	 COHM2CO2H = HCOCO2H + CO + HO2 : 	J(17) 	;")
    # Set species and new products true
    newspc[[3,4,7]] = true
    # Add species definitions
    insert!(mech,spc_end,"COHM2CO2H = IGNORE ;"); spc_end += 1; beg_rxn += 1
  end
  # Add reactions for GLYOX, if currently not in the mechanism
  if !newspc[4] && any(x->x==4,chk)
    rxn += 1
    push!(mech, "{$rxn.} 	 GLYOX = CO + CO + H2 : 	J(31) 	;")
    rxn += 1
    push!(mech, "{$rxn.} 	 GLYOX = CO + CO + HO2 + HO2 : 	J(33) 	;")
    rxn += 1
    push!(mech, "{$rxn.} 	 GLYOX = HCHO + CO : 	J(32) 	;")
    rxn += 1
    push!(mech, "{$rxn.} 	 NO3 + GLYOX = HCOCO + HNO3 : 	KNO3AL 	;")
    rxn += 1
    push!(mech, "{$rxn.} 	 OH + GLYOX = HCOCO : 	3.1D-12*EXP(340/TEMP) 	;")
    # Set species and new products true
    newspc[4:5] = true
    # Add species definitions
    insert!(mech,spc_end,"GLYOX = IGNORE ;"); spc_end += 1; beg_rxn += 1
  end

  #––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#

  # Common products from correction 1 and 2:

  # Add reactions for HCOCO, if currently not in the mechanism
  if !newspc[5] && any(x->x==5,chk)
    rxn += 1
    push!(mech, "{$rxn.} 	 HCOCO = CO + CO + HO2 : 	7.00D11*EXP(-3160/TEMP)+5.00D-12*O2 	;")
    rxn += 1
    push!(mech, "{$rxn.} 	 HCOCO = CO + OH : 	5.00D-12*O2*3.2*(1-EXP(-550/TEMP)) 	;")
    rxn += 1
    push!(mech, "{$rxn.} 	 HCOCO = HCOCO3 : 	5.00D-12*O2*3.2*EXP(-550/TEMP) 	;")
    # Set species and new products true
    newspc[5:6] = true
    # Add species definitions
    insert!(mech,spc_end,"HCOCO = IGNORE ;"); spc_end += 1; beg_rxn += 1
  end
  # Add reactions for HCOCO3, if currently not in the mechanism
  if !newspc[6] && any(x->x==6,chk)
    rxn += 1
    push!(mech, "{$rxn.} 	 HCOCO3 + HO2 = HCOCO2H + O3 : 	KAPHO2*0.15 	;")
    rxn += 1
    push!(mech, "{$rxn.} 	 HCOCO3 + HO2 = HCOCO3H : 	KAPHO2*0.41 	;")
    rxn += 1
    push!(mech, "{$rxn.} 	 HCOCO3 + HO2 = HO2 + CO + OH : 	KAPHO2*0.44 	;")
    rxn += 1
    push!(mech, "{$rxn.} 	 HCOCO3 + NO = HO2 + CO + NO2 : 	KAPNO 	;")
    rxn += 1
    push!(mech, "{$rxn.} 	 HCOCO3 + NO2 = HO2 + CO + NO3 : 	KFPAN 	;")
    rxn += 1
    push!(mech, "{$rxn.} 	 HCOCO3 + NO3 = HO2 + CO + NO2 : 	KRO2NO3*1.74 	;")
    rxn += 1
    push!(mech, "{$rxn.} 	 HCOCO3 = CO + HO2 : 	1.00D-11*0.7*RO2 	;")
    rxn += 1
    push!(mech, "{$rxn.} 	 HCOCO3 = HCOCO2H : 	1.00D-11*0.3*RO2 	;")
    # Set species and new products true
    newspc[6:8] = true
    # Add species definitions
    insert!(mech,spc_end,"HCOCO3 = IGNORE ;"); spc_end += 1; beg_rxn += 1
  end
  # Add reactions for HCOCO2H, if currently not in the mechanism
  if !newspc[7] && any(x->x==7,chk)
    rxn += 1
    push!(mech, "{$rxn.} 	 HCOCO2H = HO2 + HO2 + CO : 	J(34) 	;")
    rxn += 1
    push!(mech, "{$rxn.} 	 OH + HCOCO2H = CO + HO2 : 	1.23D-11 	;")
    # Set species true
    newspc[7] = true
    # Add species definitions
    insert!(mech,spc_end,"HCOCO2H = IGNORE ;"); spc_end += 1; beg_rxn += 1
  end
  # Add reactions for HCOCO3H, if currently not in the mechanism
  if !newspc[8] && any(x->x==8,chk)
    rxn += 1
    push!(mech, "{$rxn.} 	 HCOCO3H = HO2 + CO + OH : 	J(41)+J(15) 	;")
    rxn += 1
    push!(mech, "{$rxn.} 	 OH + HCOCO3H = HCOCO3 : 	1.58D-11 	;")
    # Set species true
    newspc[8] = true
    # Add species definitions
    insert!(mech,spc_end,"HCOCO3H = IGNORE ;"); spc_end += 1; beg_rxn += 1
  end
  # Add reactions for CH3CO3, if currently not in the mechanism
  if !newspc[9] && any(x->x==9,chk)
    rxn += 1
    push!(mech, "{$rxn.} 	 ")
    rxn += 1
    push!(mech, "{$rxn.} 	 ")
    rxn += 1
    push!(mech, "{$rxn.} 	 ")
    rxn += 1
    push!(mech, "{$rxn.} 	 ")
    rxn += 1
    push!(mech, "{$rxn.} 	 ")
    rxn += 1
    push!(mech, "{$rxn.} 	 ")
    rxn += 1
    push!(mech, "{$rxn.} 	 ")
    rxn += 1
    push!(mech, "{$rxn.} 	 ")
    # Set species and new products true
    newspc[[3,4,7]] = true
    # Add species definitions
    insert!(mech,spc_end,"COHM2CO2H = IGNORE ;"); spc_end += 1; beg_rxn += 1
  end

  #––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#

  return mech,newspc,spc_end,beg_rxn
end #function add_prod

end #module chk_mech
