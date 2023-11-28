#!/bin/bash


sed -i.bak 's/\(u_{[0-3]}\){\\left(x_{[0-3]} \\right)}/\1/g' ricci_tensor.tex
sed -i.bak 's/\\left(\\frac{d}{d x_{[0-3]}} u_{\([0-3]\)}\\right)/\\dot{u}_{\1}/g' ricci_tensor.tex
sed -i.bak 's/\\frac{d^{2}}{d x_{[0-3]}^{2}} u_{\([0-3]\)}/\\ddot{u}_{\1}/g' ricci_tensor.tex
sed -i.bak 's/\\frac{d}{d x_{[0-3]}} u_{\([0-3]\)}/\\dot{u}_{\1}/g' ricci_tensor.tex
sed -i.bak 's/- 1\.0 /-/g' ricci_tensor.tex
sed -i.bak 's/- 2 /-2/g' ricci_tensor.tex


exit 0
