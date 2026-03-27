First of all, I have to apologize for my poor habits in programming. Due to the update of the goal, the codes that were originally summarized once are summarized again. And in the codes, there are too many 'if' and few notes. That is why I write this manual.

The code files work in MATLAB and are separated in three kinds: Total running code(file name with 'ALL', usually in the root directory. May be it should be called 'main main code'), main code(it was main code before. With more goals it becomes a subroutine) and subroutine. Some subroutine like 'CG-EE-Grid-Scan' can be ran by itself.

In order to display the code in the suitable sequence, I name the function code with 'Z' to ensure them at the end, and name other subroutine with 'A' to 'H' to ensure the sequence follow the workflow. In addition, if a file is used only for CG-FD or SG-FD, the name starts from 'CG' or 'SG'. A special name, 'EE' means this file is related to 'E'. The stability checking of 'CG\_E' is based on the results from 'CG\_EE'. You need to run 'CG\_EE' with '$gammarange=0.1:0.1:0.9$' for the first time to save the workspace.

If I have time, I will upload a menual for more details.
