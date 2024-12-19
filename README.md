```bash
cmake ..  -DINCLUDE_TESTS=ON -DINCLUDE_TRAIN=ON \
        -DLPSOLVE_INCLUDE_DIR=/opt/homebrew/Cellar/lp_solve/5.5.2.11/include \
        -DLPSOLVE_LIBRARY_DIR=/opt/homebrew/Cellar/lp_solve/5.5.2.11/lib \
        -DRDKIT_INCLUDE_DIR=/opt/homebrew/Cellar/rdkit/2024.09.3_1/include/rdkit/ \
        -DRDKIT_INCLUDE_EXT_DIR=/opt/homebrew/Cellar/rdkit/2024.09.3_1/include/rdkit/GraphMol \
        -DRDKIT_LIBRARY_DIR=/opt/homebrew/Cellar/rdkit/2024.09.3_1/lib \
        -DCMAKE_CXX_STANDARD=20 \
        -DCMAKE_BUILD_TYPE=Release
```


problem is at line 795 of FragmentTreeNode.cpp where the calculated valence of the carbon is greater than allowed.
