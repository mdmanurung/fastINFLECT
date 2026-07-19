# Plot the two halves of the curve with touchlines and display angle

Internal function. Displays basic plot of curve and touchlines created
by `leastError`. Of greater use when refining Inflection Point

## Usage

``` r
angleplot(part1, part2, test1, test2, main, angle)
```

## Arguments

  - part1:
    
    Datapoints from start to kneepoint

  - part2:
    
    Datapoints from kneepoint to end

  - test1:
    
    Fitted line over part1

  - test2:
    
    Fitted line over part2

  - main:
    
    Title of the plot

  - angle:
    
    Angle between test1 and test2, obtained through `LinesAngles`
    function from `LearnGeom`

## Value

None, plot is generated using the basic R plot function.

## See also

`INFLECT` , `QC.to.curve`,`leastError`
