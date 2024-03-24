#### Introduction: 

This is repo is a test repo for simulating basic model of discrete elastic rod. It is a modifictaion  to polyscope where I change some of UI.For now you can simply create a test Bezier curve and generate tangent, normal, binormal and bishop frame on each edges. 

#### Download & Compilation:

Add –recursive when git cloning the repo. 

On mac you can compile like 

```
cd exmaples/demo-app
mkdir build 
cd build 
cmake ..
make -j8
```

Thanks to polyscope. Few dependencies is required to compile this project. And we didnt involve halfedge data structure for now. 

The compilation on windows is not tested. 

##### Future work

This is a simple test and we will migrate the main work to other framwork. 