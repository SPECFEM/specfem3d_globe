# GLAD-M25 model

You can download all binary files from the data repository [SPECFEM specfem-data](https://github.com/SPECFEM/specfem-data)<br>
This is due to the file model sizes, that we provide them in a separate data repository.
Please check out the data repository for further infos.

As an alternative, you can also download the full original repository [gladm25](https://github.com/caiociardelli/gladm25)
into this folder and setup the model by the following steps:

1. clone repository into this folder:
```
git clone https://github.com/caiociardelli/gladm25
```

2. create the two symbolic links for crust & mantle file needed for the code to find all corresponding model files:
```
ln -s gladm25/crust/
ln -s gladm25/mantle/
```

3. extract the tar files in the crust folder:
```
for file in crust/*tar.xz; do tar -xvJf "$file" -C crust/; done
```

You're all done. Now, in Par_file you can use `MODEL = glad_bkmns`.

#### gladm15
Note that in case you would be interested in using model `gladm15` rather than `gladm25`, you repeat the above steps and just clone from repository [gladm15](https://github.com/caiociardelli/gladm15):
```
git clone https://github.com/caiociardelli/gladm15
```
and correspondingly
```
ln -s gladm15/crust/
ln -s gladm15/mantle/
```
to replace the `gladm25` model files. Everything else can stay the same.
