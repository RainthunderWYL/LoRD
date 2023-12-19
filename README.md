# LoRD

**LoRD**, short for **Lo**cate **R**econnection **D**istribution, is a toolkit designed to identify the locations and structures of 3D magnetic reconnection within discrete magnetic field data.

**LoRD** contains three main functions (for now), namely, **ARD** (Analyze Reconnection Distribution), **ANP** (Analyze Null Points), and **APNP** (Analyze Projected Null Points).

**ARD** locates the grids undergoing reconnection without null points and also recognizes the local configurations of reconnection sites.

**ANP** locates and classifies the 3D null points based on the theory by Parnel et al. 1996.

**APNP** analyzes the 2D neutral points projected on a plane near a cell which can be defined by users (please don't use this function unless you know its details!).

**Please refer to the [Wiki](https://github.com/RainthunderWYL/LoRD/wiki) page for details.**
