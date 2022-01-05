These files are an implementation of the algorithm described in Fraser and Schwartz, 
"Recording from the same neurons chronically in motor cortex", Journal of Neurophysiology 
2012 (in press).

The main function is "unitIdentification".  The main output is "survival", a cell array of
matrices describing which recorded units represent neurons that survive from one recording 
session to the next.  

The function "identifiers" converts "survival" into a different form which may be more
convenient in some applications.  Instead of returning a matrix of 1s and 0s indicating
survival across sessions, "identifiers" assigns a unique identifier to every neuron that
is tracked for one or more sessions, and then labels every recorded unit with this global
identifier.

George Fraser
12/9/2011
fraser.george.w@gmail.com