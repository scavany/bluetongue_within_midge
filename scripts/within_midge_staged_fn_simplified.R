## within midge model just virus
wv.BTV = function(t, state, parameters)
{
  with(as.list(c(state, parameters)),{
    dV.m = -c.m * V.m
    dV.h = p.m + k * p.s - c.s * V.h 
    list(c(dV.m,dV.h))
  })
}

## within midge model intrathoracic (no eclipse)
wv.BTV.intrathoracic = function(t, state, parameters)
{
  with(as.list(c(state, parameters)),{
    dV.h = p.s - c.s * V.h # all infected cells produce virus
    list(c(dV.h))
  })
}
