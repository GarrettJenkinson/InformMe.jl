# This file is a modified version of unique.m from the octave language that was downloaded here:
# http://code.metager.de/source/xref/gnu/octave/scripts/set/unique.m
# It has been modified by Garrett Jenkinson to work in Julia.
# It came with the following header, and thus is distributed under GPLv3.

## Copyright (C) 2000-2013 Paul Kienzle
## Copyright (C) 2008-2009 Jaroslav Hajek
##
## This file is part of Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or (at
## your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {Function File} {} unique (@var{x})
## @deftypefnx {Function File} {} unique (@var{x}, "rows")
## @deftypefnx {Function File} {} unique (@dots{}, "first")
## @deftypefnx {Function File} {} unique (@dots{}, "last")
## @deftypefnx {Function File} {[@var{y}, @var{i}, @var{j}] =} unique (@dots{})
## Return the unique elements of @var{x}, sorted in ascending order.
## If the input @var{x} is a vector then the output is also a vector with the
## same orientation (row or column) as the input.  For a matrix input the
## output is always a column vector.  @var{x} may also be a cell array of
## strings.
##
## If the optional argument @qcode{"rows"} is supplied, return the unique
## rows of @var{x}, sorted in ascending order.
##
## If requested, return index vectors @var{i} and @var{j} such that
## @code{x(i)==y} and @code{y(j)==x}.
##
## Additionally, if @var{i} is a requested output then one of @qcode{"first"} or
## @qcode{"last"} may be given as an input.  If @qcode{"last"} is specified,
## return the highest possible indices in @var{i}, otherwise, if @qcode{"first"}
## is specified, return the lowest.  The default is @qcode{"last"}.
## @seealso{union, intersect, setdiff, setxor, ismember}
## @end deftypefn

function uniqueOct(x)
  n = length(x)

  y = deepcopy(x)
  ## Special cases 0 and 1
  if (n == 0)
    i = j = []
    return y, i, j
  elseif (n == 1)
    i = j = 1
    return y, i, j
  end

  i=sortperm(y)
  sort!(y)

  match = falses(n-1)
  intMatch = zeros(Int,n-1)
  @inbounds for index=1:(n-1)
    match[index] = y[index] == y[index+1]
    if match[index]
      intMatch[index]=1
    end
  end
  idx = findall(match)

  j = deepcopy(i)

  j[i] = cumsum(cat(1,(1 .- intMatch),dims=1))

  @inbounds for index=length(idx):-1:1
    splice!(y,idx[index])
    splice!(i,idx[index]+1)
  end

  return y, i, j

end
