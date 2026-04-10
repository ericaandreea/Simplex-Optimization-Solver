A simple command-line application that solves linear programming problems using the Simplex Algorithm.

Technical Features:

- Uses the Simplex method with the Big-M technique to handle constraints of type =, >=, and <=.
- Displays the full simplex tableau at each iteration, providing clear traceability of variables entering and leaving the basis.
- Includes a post-calculation verification module to ensure the solution satisfies all initial constraints.
   
![rezultat](https://github.com/user-attachments/assets/3bcbba48-b9c5-47f1-8d80-63b113006b6f)

How to use:

- Compile the file: g++ main.cpp -o solver
- Run the executable: ./solver
- Enter the problem data according to the instructions displayed in the console.
