#! /usr/bin/env python
from PyTrilinos import Epetra, Triutils

Epetra.Init()
Comm = Epetra.PyComm()
n = 1000000 * Comm.NumProc()
Gallery = Triutils.CrsMatrixGallery("laplace_1d", Comm)
Gallery.Set("problem_size", n)
Matrix = Gallery.GetMatrix()
LHS = Gallery.GetStartingSolution()
RHS = Gallery.GetRHS()
Time = Epetra.Time(Comm)

for i in xrange(100):
  Matrix.Multiply(False, LHS, RHS)

print Time.ElapsedTime()
Epetra.Finalize()
