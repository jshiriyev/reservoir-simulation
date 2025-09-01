# porsim  

**porsim** is a Python package for **reservoir simulation workflows**, designed with education and research in mind. It implements core simulation concepts and exercises following Balhoff’s *Reservoir Simulation* book, and provides a foundation for benchmarking and experimentation in porous media flow.  

The package is continuously developed and improved every semester as part of a graduate/undergraduate reservoir simulation course, making it both a **learning resource** and a **research toolkit**.  

---

## ✨ Features  

- **Tutorial Workflows** based on Balhoff’s book  
- **Single-phase flow**: one-dimensional analytical benchmarks  
- **Two-phase immiscible flow**: Buckley–Leverett solutions  
- **Non-linear solvers**: Picard Iteration methods for convergence  
- **Benchmarking utilities** for comparing analytical vs. numerical results  
- **Educational focus**: well-documented code and examples for classroom use  
- **Work in progress**: continuously expanded with new features each semester  

---

## 📚 Examples  

- **Single-phase flow benchmark**  
  Generate 1D analytical solutions to validate numerical schemes.  

- **Buckley–Leverett problem**  
  Explore immiscible displacement in porous media with fractional flow analysis.  

- **Non-linear solvers**  
  Apply Picard iteration to handle non-linearities in flow equations.  

```python
import porsim  

# Example: run a single-phase benchmark  
sim = porsim.Solver(length=100, nx=50)  
sim.run()
sim.plot_solution()  
```

---

## 🚧 Development Status  

This package is **actively developed** as part of a reservoir simulation course.  
New methods and workflows are added each semester, making it a living repository of practical examples.  

Planned future extensions:  
- Multidimensional flow  
- Multiphase capillary pressure & relative permeability  
- Advanced solvers and preconditioners  
- Coupled well–reservoir models  

---

## 🔧 Installation  

```bash
pip install porsim
```

Or from source:  

https://github.com/jshiriyev/reservoir-simulation.git

```bash
git clone https://github.com/jshiriyev/reservoir-simulation.git
cd reservoir-simulation
pip install -e .
```

---

## 🤝 Contributing  

Contributions are welcome, especially improvements to examples, documentation, and solver methods.  
This project evolves through teaching and collaboration — feedback from students, researchers, and practitioners is highly valued.  

---

## 📬 Contact  

Maintained by **Javid Shiriyev**  
- 📧 [shiriyevcavid@gmail.com](mailto:shiriyevcavid@gmail.com)  
- 🔗 [LinkedIn](https://www.linkedin.com/in/jshiriyev/)  

