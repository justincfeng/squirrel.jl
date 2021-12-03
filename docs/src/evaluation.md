# Evaluation

The functions described here are used to evaluate the locator functions
for the `squirrel.jl` code.

## Structs

```@docs
squirrel.seval.TestCases
```

```@docs
squirrel.seval.TestData
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

## Evaluation

The following is the main evaluation function. It returns a quantity of
datatype TestData.

```@docs
squirrel.seval.main
```

