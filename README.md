```bash
cmake ..  -DINCLUDE_TESTS=ON -DINCLUDE_TRAIN=ON \                                                         28s  base  !4/7
        -DLPSOLVE_INCLUDE_DIR=/opt/homebrew/Cellar/lp_solve/5.5.2.11/include \
        -DLPSOLVE_LIBRARY_DIR=/opt/homebrew/Cellar/lp_solve/5.5.2.11/lib \
        -DRDKIT_INCLUDE_DIR=/opt/homebrew/Cellar/rdkit/2024.09.3_1/include/rdkit/ \
        -DRDKIT_INCLUDE_EXT_DIR=/opt/homebrew/Cellar/rdkit/2024.09.3_1/include/rdkit/GraphMol \
        -DRDKIT_LIBRARY_DIR=/opt/homebrew/Cellar/rdkit/2024.09.3_1/lib \
        -DCMAKE_CXX_STANDARD=17 \
        -DCMAKE_BUILD_TYPE=Release
```
