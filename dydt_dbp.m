## Copyright (C) 2023 Arghya
## 
## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see
## <https://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {} {@var{retval} =} dydt (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Arghya <arghya@arghya-Inspiron-5558>
## Created: 2023-04-15

function stvd = dydt_dbp (t,stv)
  m1=1;m2=1;L1=0.5;L2=0.5;g=10;
  q1=stv(1);q2=stv(2);u1=stv(3);u2=stv(4);
  mm = [m1*L1^2+m2*L1^2 m2*L1*L2*cos(q1-q2); m2*L1*L2*cos(q1-q2) m2*L2^2];
  forcing = [-L1*(L2*m2*sin(q1-q2)*u2^2+g*m1*sin(q1)+g*m2*sin(q1));...
                                 L2*m2*(L1*sin(q1-q2)*u1^2-g*sin(q2))];
  amat=[eye(2) zeros(2); zeros(2) mm];
  bvec=[[u1;u2]; forcing];
  stvd=amat\bvec;

endfunction
