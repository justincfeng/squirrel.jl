# Evaluation functions

The functions described here are contained in the `seval` submodule.
These functions are encapsulated so that they provide a somewhat
independent evaluation of the locator functions in the `squirrel.jl`
code. At present, the parameters of the evaluation functions have values
corresponding to terrestrial positioning.

## Basic strategy

The functions in `seval`, particularly `seval.gen` and `seval.main`, may
be used to evaluate the accuracy of the locator functions implemented in
`squirrel.jl` according to the following strategy:

1.  Generate (stochastically) a target point `Xtar`.

2.  Generate `ne` emission points `X` on the past light cone of `Xtar`
    by solving the geodesic equation for some metric `g`.

3.  Feed the emission points `X` and metric `g` into the
    `squirrel.locator` or `squirrel.locator4` function and compare with
    target point `Xtar`.

The function `seval.gen` is implements steps 1. and 2., and the function
`seval.main` implements step 3.

## Structs and format conversion functions

The evaluation submodule `seval` defines two datatypes, one for the
generation of test cases (`TestCases`), and one for storing the results
of the evaluation (`TestData`).

```@docs
squirrel.seval.TestCases
```

```@docs
squirrel.seval.TestData
```

The functions `seval.gen` and `seval.main` can generate large amounts of
data. The data may be saved to a file using Julia's built in serializer,
but it is recommended that the data be converted to a tuple. The following
functions may be used to convert the between tuples and the 
`TestCases`/`TestData` datatypes:

```@docs
squirrel.seval.tc2tup
```

```@docs
squirrel.seval.td2tup
```

```@docs
squirrel.seval.tup2tc
```

```@docs
squirrel.seval.tup2td
```

In some cases, one may wish to change to a higher precision floating point type for the generated samples. For this purpose, the following function is provided:

```@docs
squirrel.seval.tcfl
```

## Test case generation

### Point generation

Emission points are generated on the past light cone of a target point.
The idea is to stochastically generate a target point `Xtar`, and then
generate initial data for `N` past-directed null geodesics, which are
then integrated to obtain the emission points.

```@docs
squirrel.seval.pgen
```

### Multiple case generator

The following function generates `N` sets of target points `Xtar` and
emission points `X`. It returns a quantity of datatype TestCases.

```@docs
squirrel.seval.gen
```

## Evaluation function

The following is the main evaluation function. It returns a quantity of
datatype TestData.

```@docs
squirrel.seval.main
```

