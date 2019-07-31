# ph-calculator
The program to predict neutralization titration reaction which studied and developed in "Haksung High School" for the "2019 Science Exhibition" in South Korea.

## Patch Note

### v1.0.0-alpha
- Initial release

### v1.0.1-alpha
- When [OH<sup>-</sup>] remains, calculate by assuming [H<sup>+</sup>] = -[OH<sup>-</sup>])

### v1.1.0-alpha
- Can calculate weak acids, weak bases and polyprotic acids.

### v1.1.1-alpha
- Add graph-generator.

### v1.1.2-alpha
- Add equivalence-point-finder.

### v1.1.3-alpha
- Rename the equivalence-point-finder function to recipe-finder function.
- Modify the output messages.

### v1.1.4-alpha
- Solve the infinite loop problem that occurred in recipe-finder function.
- Improve the accuracy of recipe-finder function.

### v1.1.5-alpha
- Remove the danger of buffer overflow.
- Modify the output message of recipe-finder.

### v1.1.6-alpha
- Eliminate the risk of potential errors in the graph-generator function.

### v1.1.7-alpha
- Separate the functions for calculation as a library.
- Make the CalcInitialH function executed inside the CalculatePH function.
- Eliminate some risks of potential errors in the recipe-finder function.
- Change the path of the solute.pcd and result.pcd.
- A little refactoring.

### v1.1.8-alpha
- Eliminate some risks of potential errors while loading DB.

### v1.1.9-alpha
- Eliminate the risk of potential error while generating graph.

### v1.1.10-alpha
- Fix a small error with escape condition in loop while generating graph.
- Add some informational messages while calculating.
