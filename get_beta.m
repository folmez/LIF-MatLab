function b = get_beta (ge, gi, gL, vE, vI, vL, extI)

b = ge.*vE + gi.*vI + gL.*vL + extI;

end
