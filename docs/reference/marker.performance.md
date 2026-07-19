# Plot marker performance across metaclustering results

Plots the marker performance across metaclusterings in a boxplot.
Performance score is the percentage of clusters where a marker passed
the [dip.test](https://rdrr.io/pkg/diptest/man/dip.test.html) and IQR
test per metaclustering. Information is displayed in a boxplot, color
gradient denotes `set.i` .

## Usage

``` r
marker.performance(inflect.results, ggtitle = NULL, markers = NULL)
```

## Arguments

  - inflect.results:
    
    A inflect.results object resulting from `INFLECT`.

  - ggtitle:
    
    Optional. Character string to plot as title.

  - markers:
    
    Which markers should be included in the plot? A vector of strings
    matching evaluated marker names. If `NULL`, all evaluated markers
    are displayed.

## Value

`list` with 2 items. First is the melted dataframe with the marker
performance (\\ performance on the y-axis, markers on the x-axis and the
color scale denoting the amount of metaclusters evaluated.

## Examples

``` r
# Read in FlowSOM object from file. Downsampled clustering result of Levine32 dataset clustering.
# SOM-clustered to 375 clusters.
flowsom <- system.file("extdata", "Levine32sample.Rdata", package="fastINFLECT")
load(flowsom)
inflect.results<- INFLECT(FlowSOM.results= dataset, set.i= 5:12, multicore=FALSE, zeroes.in=FALSE)
#> Warning: Warning: Number of points in set.i is low, beware of noisy diagnostic curves 

# Display diagnostic graph
inflect.results$ggplot


# Now check marker performance for all markers
marker.performance(inflect.results, ggtitle= "Levine32sample", markers=NULL)
#> $marker.dataframe
#>      i           Marker Performance
#> 1    5     CD3(Er170)Di    60.00000
#> 2    6     CD3(Er170)Di    83.33333
#> 3    7     CD3(Er170)Di    85.71429
#> 4    8     CD3(Er170)Di    87.50000
#> 5    9     CD3(Er170)Di    77.77778
#> 6   10     CD3(Er170)Di    80.00000
#> 7   11     CD3(Er170)Di    81.81818
#> 8   12     CD3(Er170)Di    83.33333
#> 9    5     CD4(Nd145)Di   100.00000
#> 10   6     CD4(Nd145)Di   100.00000
#> 11   7     CD4(Nd145)Di   100.00000
#> 12   8     CD4(Nd145)Di    87.50000
#> 13   9     CD4(Nd145)Di    88.88889
#> 14  10     CD4(Nd145)Di    90.00000
#> 15  11     CD4(Nd145)Di    90.90909
#> 16  12     CD4(Nd145)Di    91.66667
#> 17   5     CD7(Dy162)Di   100.00000
#> 18   6     CD7(Dy162)Di   100.00000
#> 19   7     CD7(Dy162)Di   100.00000
#> 20   8     CD7(Dy162)Di   100.00000
#> 21   9     CD7(Dy162)Di   100.00000
#> 22  10     CD7(Dy162)Di   100.00000
#> 23  11     CD7(Dy162)Di   100.00000
#> 24  12     CD7(Dy162)Di   100.00000
#> 25   5     CD8(Nd146)Di    80.00000
#> 26   6     CD8(Nd146)Di    83.33333
#> 27   7     CD8(Nd146)Di    85.71429
#> 28   8     CD8(Nd146)Di   100.00000
#> 29   9     CD8(Nd146)Di   100.00000
#> 30  10     CD8(Nd146)Di   100.00000
#> 31  11     CD8(Nd146)Di   100.00000
#> 32  12     CD8(Nd146)Di   100.00000
#> 33   5   CD11b(Nd144)Di    80.00000
#> 34   6   CD11b(Nd144)Di    83.33333
#> 35   7   CD11b(Nd144)Di    85.71429
#> 36   8   CD11b(Nd144)Di    87.50000
#> 37   9   CD11b(Nd144)Di    88.88889
#> 38  10   CD11b(Nd144)Di    90.00000
#> 39  11   CD11b(Nd144)Di    90.90909
#> 40  12   CD11b(Nd144)Di    91.66667
#> 41   5   CD11c(Tb159)Di   100.00000
#> 42   6   CD11c(Tb159)Di   100.00000
#> 43   7   CD11c(Tb159)Di   100.00000
#> 44   8   CD11c(Tb159)Di   100.00000
#> 45   9   CD11c(Tb159)Di   100.00000
#> 46  10   CD11c(Tb159)Di   100.00000
#> 47  11   CD11c(Tb159)Di   100.00000
#> 48  12   CD11c(Tb159)Di   100.00000
#> 49   5    CD13(Er168)Di   100.00000
#> 50   6    CD13(Er168)Di   100.00000
#> 51   7    CD13(Er168)Di   100.00000
#> 52   8    CD13(Er168)Di   100.00000
#> 53   9    CD13(Er168)Di   100.00000
#> 54  10    CD13(Er168)Di   100.00000
#> 55  11    CD13(Er168)Di   100.00000
#> 56  12    CD13(Er168)Di   100.00000
#> 57   5    CD14(Gd156)Di   100.00000
#> 58   6    CD14(Gd156)Di   100.00000
#> 59   7    CD14(Gd156)Di   100.00000
#> 60   8    CD14(Gd156)Di   100.00000
#> 61   9    CD14(Gd156)Di   100.00000
#> 62  10    CD14(Gd156)Di   100.00000
#> 63  11    CD14(Gd156)Di   100.00000
#> 64  12    CD14(Gd156)Di   100.00000
#> 65   5    CD15(Dy164)Di   100.00000
#> 66   6    CD15(Dy164)Di   100.00000
#> 67   7    CD15(Dy164)Di   100.00000
#> 68   8    CD15(Dy164)Di   100.00000
#> 69   9    CD15(Dy164)Di   100.00000
#> 70  10    CD15(Dy164)Di   100.00000
#> 71  11    CD15(Dy164)Di   100.00000
#> 72  12    CD15(Dy164)Di   100.00000
#> 73   5    CD16(Ho165)Di   100.00000
#> 74   6    CD16(Ho165)Di    83.33333
#> 75   7    CD16(Ho165)Di    85.71429
#> 76   8    CD16(Ho165)Di    87.50000
#> 77   9    CD16(Ho165)Di    88.88889
#> 78  10    CD16(Ho165)Di    90.00000
#> 79  11    CD16(Ho165)Di    90.90909
#> 80  12    CD16(Ho165)Di    91.66667
#> 81   5    CD19(Nd142)Di    80.00000
#> 82   6    CD19(Nd142)Di    83.33333
#> 83   7    CD19(Nd142)Di    85.71429
#> 84   8    CD19(Nd142)Di    87.50000
#> 85   9    CD19(Nd142)Di    88.88889
#> 86  10    CD19(Nd142)Di    90.00000
#> 87  11    CD19(Nd142)Di   100.00000
#> 88  12    CD19(Nd142)Di   100.00000
#> 89   5    CD20(Sm147)Di   100.00000
#> 90   6    CD20(Sm147)Di   100.00000
#> 91   7    CD20(Sm147)Di   100.00000
#> 92   8    CD20(Sm147)Di   100.00000
#> 93   9    CD20(Sm147)Di   100.00000
#> 94  10    CD20(Sm147)Di   100.00000
#> 95  11    CD20(Sm147)Di   100.00000
#> 96  12    CD20(Sm147)Di   100.00000
#> 97   5    CD22(Nd143)Di   100.00000
#> 98   6    CD22(Nd143)Di   100.00000
#> 99   7    CD22(Nd143)Di   100.00000
#> 100  8    CD22(Nd143)Di   100.00000
#> 101  9    CD22(Nd143)Di   100.00000
#> 102 10    CD22(Nd143)Di   100.00000
#> 103 11    CD22(Nd143)Di   100.00000
#> 104 12    CD22(Nd143)Di   100.00000
#> 105  5    CD33(Gd158)Di   100.00000
#> 106  6    CD33(Gd158)Di   100.00000
#> 107  7    CD33(Gd158)Di   100.00000
#> 108  8    CD33(Gd158)Di   100.00000
#> 109  9    CD33(Gd158)Di   100.00000
#> 110 10    CD33(Gd158)Di   100.00000
#> 111 11    CD33(Gd158)Di   100.00000
#> 112 12    CD33(Gd158)Di   100.00000
#> 113  5    CD34(Nd148)Di    80.00000
#> 114  6    CD34(Nd148)Di    83.33333
#> 115  7    CD34(Nd148)Di   100.00000
#> 116  8    CD34(Nd148)Di   100.00000
#> 117  9    CD34(Nd148)Di   100.00000
#> 118 10    CD34(Nd148)Di   100.00000
#> 119 11    CD34(Nd148)Di   100.00000
#> 120 12    CD34(Nd148)Di   100.00000
#> 121  5    CD38(Er167)Di    80.00000
#> 122  6    CD38(Er167)Di    66.66667
#> 123  7    CD38(Er167)Di    85.71429
#> 124  8    CD38(Er167)Di    87.50000
#> 125  9    CD38(Er167)Di    88.88889
#> 126 10    CD38(Er167)Di    90.00000
#> 127 11    CD38(Er167)Di    90.90909
#> 128 12    CD38(Er167)Di    91.66667
#> 129  5    CD41(Lu175)Di   100.00000
#> 130  6    CD41(Lu175)Di   100.00000
#> 131  7    CD41(Lu175)Di   100.00000
#> 132  8    CD41(Lu175)Di   100.00000
#> 133  9    CD41(Lu175)Di   100.00000
#> 134 10    CD41(Lu175)Di   100.00000
#> 135 11    CD41(Lu175)Di   100.00000
#> 136 12    CD41(Lu175)Di   100.00000
#> 137  5    CD44(Er166)Di   100.00000
#> 138  6    CD44(Er166)Di   100.00000
#> 139  7    CD44(Er166)Di   100.00000
#> 140  8    CD44(Er166)Di   100.00000
#> 141  9    CD44(Er166)Di   100.00000
#> 142 10    CD44(Er166)Di   100.00000
#> 143 11    CD44(Er166)Di   100.00000
#> 144 12    CD44(Er166)Di   100.00000
#> 145  5    CD45(Sm154)Di   100.00000
#> 146  6    CD45(Sm154)Di   100.00000
#> 147  7    CD45(Sm154)Di   100.00000
#> 148  8    CD45(Sm154)Di   100.00000
#> 149  9    CD45(Sm154)Di   100.00000
#> 150 10    CD45(Sm154)Di   100.00000
#> 151 11    CD45(Sm154)Di   100.00000
#> 152 12    CD45(Sm154)Di   100.00000
#> 153  5  CD45RA(La139)Di   100.00000
#> 154  6  CD45RA(La139)Di   100.00000
#> 155  7  CD45RA(La139)Di   100.00000
#> 156  8  CD45RA(La139)Di   100.00000
#> 157  9  CD45RA(La139)Di   100.00000
#> 158 10  CD45RA(La139)Di   100.00000
#> 159 11  CD45RA(La139)Di   100.00000
#> 160 12  CD45RA(La139)Di   100.00000
#> 161  5    CD47(Gd160)Di   100.00000
#> 162  6    CD47(Gd160)Di   100.00000
#> 163  7    CD47(Gd160)Di   100.00000
#> 164  8    CD47(Gd160)Di   100.00000
#> 165  9    CD47(Gd160)Di   100.00000
#> 166 10    CD47(Gd160)Di   100.00000
#> 167 11    CD47(Gd160)Di   100.00000
#> 168 12    CD47(Gd160)Di   100.00000
#> 169  5   CD49d(Yb172)Di   100.00000
#> 170  6   CD49d(Yb172)Di   100.00000
#> 171  7   CD49d(Yb172)Di   100.00000
#> 172  8   CD49d(Yb172)Di   100.00000
#> 173  9   CD49d(Yb172)Di   100.00000
#> 174 10   CD49d(Yb172)Di   100.00000
#> 175 11   CD49d(Yb172)Di   100.00000
#> 176 12   CD49d(Yb172)Di   100.00000
#> 177  5    CD61(Tm169)Di   100.00000
#> 178  6    CD61(Tm169)Di   100.00000
#> 179  7    CD61(Tm169)Di   100.00000
#> 180  8    CD61(Tm169)Di   100.00000
#> 181  9    CD61(Tm169)Di   100.00000
#> 182 10    CD61(Tm169)Di   100.00000
#> 183 11    CD61(Tm169)Di   100.00000
#> 184 12    CD61(Tm169)Di   100.00000
#> 185  5    CD64(Yb176)Di   100.00000
#> 186  6    CD64(Yb176)Di   100.00000
#> 187  7    CD64(Yb176)Di   100.00000
#> 188  8    CD64(Yb176)Di   100.00000
#> 189  9    CD64(Yb176)Di   100.00000
#> 190 10    CD64(Yb176)Di    90.00000
#> 191 11    CD64(Yb176)Di    90.90909
#> 192 12    CD64(Yb176)Di    91.66667
#> 193  5   CD117(Yb171)Di   100.00000
#> 194  6   CD117(Yb171)Di   100.00000
#> 195  7   CD117(Yb171)Di   100.00000
#> 196  8   CD117(Yb171)Di   100.00000
#> 197  9   CD117(Yb171)Di   100.00000
#> 198 10   CD117(Yb171)Di   100.00000
#> 199 11   CD117(Yb171)Di   100.00000
#> 200 12   CD117(Yb171)Di   100.00000
#> 201  5   CD123(Eu151)Di    60.00000
#> 202  6   CD123(Eu151)Di    66.66667
#> 203  7   CD123(Eu151)Di    71.42857
#> 204  8   CD123(Eu151)Di    75.00000
#> 205  9   CD123(Eu151)Di    77.77778
#> 206 10   CD123(Eu151)Di    80.00000
#> 207 11   CD123(Eu151)Di    90.90909
#> 208 12   CD123(Eu151)Di    91.66667
#> 209  5   CD133(Pr141)Di   100.00000
#> 210  6   CD133(Pr141)Di   100.00000
#> 211  7   CD133(Pr141)Di   100.00000
#> 212  8   CD133(Pr141)Di   100.00000
#> 213  9   CD133(Pr141)Di   100.00000
#> 214 10   CD133(Pr141)Di   100.00000
#> 215 11   CD133(Pr141)Di   100.00000
#> 216 12   CD133(Pr141)Di   100.00000
#> 217  5 CD235ab(Sm152)Di   100.00000
#> 218  6 CD235ab(Sm152)Di   100.00000
#> 219  7 CD235ab(Sm152)Di   100.00000
#> 220  8 CD235ab(Sm152)Di   100.00000
#> 221  9 CD235ab(Sm152)Di   100.00000
#> 222 10 CD235ab(Sm152)Di   100.00000
#> 223 11 CD235ab(Sm152)Di   100.00000
#> 224 12 CD235ab(Sm152)Di   100.00000
#> 225  5   CD321(Eu153)Di   100.00000
#> 226  6   CD321(Eu153)Di   100.00000
#> 227  7   CD321(Eu153)Di   100.00000
#> 228  8   CD321(Eu153)Di   100.00000
#> 229  9   CD321(Eu153)Di   100.00000
#> 230 10   CD321(Eu153)Di   100.00000
#> 231 11   CD321(Eu153)Di   100.00000
#> 232 12   CD321(Eu153)Di   100.00000
#> 233  5   CXCR4(Sm149)Di   100.00000
#> 234  6   CXCR4(Sm149)Di   100.00000
#> 235  7   CXCR4(Sm149)Di   100.00000
#> 236  8   CXCR4(Sm149)Di   100.00000
#> 237  9   CXCR4(Sm149)Di   100.00000
#> 238 10   CXCR4(Sm149)Di   100.00000
#> 239 11   CXCR4(Sm149)Di   100.00000
#> 240 12   CXCR4(Sm149)Di   100.00000
#> 241  5    Flt3(Nd150)Di   100.00000
#> 242  6    Flt3(Nd150)Di   100.00000
#> 243  7    Flt3(Nd150)Di   100.00000
#> 244  8    Flt3(Nd150)Di   100.00000
#> 245  9    Flt3(Nd150)Di    88.88889
#> 246 10    Flt3(Nd150)Di    90.00000
#> 247 11    Flt3(Nd150)Di    90.90909
#> 248 12    Flt3(Nd150)Di    91.66667
#> 249  5  HLA-DR(Yb174)Di    40.00000
#> 250  6  HLA-DR(Yb174)Di    50.00000
#> 251  7  HLA-DR(Yb174)Di    42.85714
#> 252  8  HLA-DR(Yb174)Di    50.00000
#> 253  9  HLA-DR(Yb174)Di    44.44444
#> 254 10  HLA-DR(Yb174)Di    50.00000
#> 255 11  HLA-DR(Yb174)Di    54.54545
#> 256 12  HLA-DR(Yb174)Di    50.00000
#> 
#> $plot
#> Warning: The following aesthetics were dropped during statistical transformation:
#> colour.
#> ℹ This can happen when ggplot fails to infer the correct grouping structure in
#>   the data.
#> ℹ Did you forget to specify a `group` aesthetic or to convert a numerical
#>   variable into a factor?

#> 

```
