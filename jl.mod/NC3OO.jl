"""
# Module NC3OO

Correct MCMv3.3.1 mechanism for NC3OO species.

# Functions (public)
- del_spc
- replNC3OO
"""
module NC3OO

export del_spc, replNC3OO


"""
    del_spc(mech,spc_end,beg_rxn)

Delete species definition of NC3OO in `mech` and adjust parameter `spc_end`
and `beg_rxn` with index of last species definition and first reaction, respectively.

Returns mech, spc_end, and beg_rxn.
"""
function del_spc(mech,spc_end,beg_rxn)

  for i = length(mech):-1:1
    if ismatch(r"NC3OO.*=",mech[i])
      if ismatch(r"NC3OO.*IGNORE",mech[i])  beg_rxn -= 1; spc_end -= 1  end
      deleteat!(mech,i)
    end
  end

  return mech,spc_end,beg_rxn

end #function del_spc


"""
    replNC3OO(mech)

Replace `NC3OO` with `NOAOO` in mechanism `mech`.


Returns mech.
"""
function replNC3OO(mech)
  for i in 1:length(mech)
    mech[i] = replace(mech[i],"NC3OO","NOAOO")
  end

  return mech
end #function replNC3OO

end #module NC3OO
