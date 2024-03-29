PARTITION_X Creates sub-partitions of response matrix for quadratic
 extrapolation bias correction.

   ------
   SYNTAX
   ------
   PARTITION_X can be used in two different ways
 
       Xpart = partition_X(X, nt, Npart, partIdx, ntPart);
 
   or
 
       [Xpart, ntPart, totNtPart] = partition_X(X, nt, Npart, partIdx);

   In the first case it is assumed that the number of trial per stimulus,
   ntPart, of the final partition matrix is known. In the second case it
   is assumed that ntPart is not know and needs to be computed by the
   function itself (see also the examples below). This dual behavior can
   be used to save time when different partitions corresponding to the 
   same Npart value need to be computed.

   ---------
   ARGUMENTS
   ---------
   X       - response matrix
   nt      - number of trials per stimulus
   Npart   - number of partitions
   partIdx - index of the partition to be extracted
   ntPart  - number of trials per stimulus of the partition matrix. In 
             Matlab's language it would be ntPart = floor(nt./Npart)

   ------
   OUTPUT
   ------
   Xpart     - partition matrix
   npart     - number of trials per stimulus of the partition matrix. In 
               Matlab's language it would be ntPart = floor(nt./Npart)
   totNtPart - total number of trials of the partition matrix. In Matlab's
               language it would be totNtPart = sum(ntPart)

   --------
   EXAMPLES
   --------
   The command

       [Xpart, ntPart, totNtPart] = partition_X(X, nt, 4, 2);

   provides the second of the four possible partitions of X together with
   the number of trials per stimulus and the total number of trials of the
   partition matrix. Alternatively

       Xpart = partition_X(X, nt, 4, 2, ntPart);
  
   also provides the second of the four possible partitions of X. In this
   case, however, the array ntPart is provided as an input to speed up
   computations.

   Copyright (C) 2010 Cesare Magri
   Version: 1.0.2

 -------
 LICENSE
 -------
 This software is distributed free under the condition that:

 1. it shall not be incorporated in software that is subsequently sold;

 2. the authorship of the software shall be acknowledged and the following
    article shall be properly cited in any publication that uses results
    generated by the software:

      Magri C, Whittingstall K, Singh V, Logothetis NK, Panzeri S: A
      toolbox for the fast information analysis of multiple-site LFP, EEG
      and spike train recordings. BMC Neuroscience 2009 10(1):81;

 3.  this notice shall remain in place in each source file.

 ----------
 DISCLAIMER
 ----------
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.