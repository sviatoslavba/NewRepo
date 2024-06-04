using System;
using System.IO;
using Microsoft.Win32;
using System.Windows;
using System.Windows.Controls;
using System.Reflection.Emit;

namespace WpfC1
{
    public partial class CalculatorWindow : Window
    {
        private CalculatorManager calculatorManager;

        public CalculatorWindow()
        {
            InitializeComponent();
            calculatorManager = new CalculatorManager();
        }

        private void ButtonSystemPow_Click(object sender, RoutedEventArgs e)
        {
            calculatorManager.SelectedSystem = new PowSystemEquations();
            TextBoxL.IsEnabled = true;
            TextBoxL.Text = "";
        }

        private void ButtonSystemCos_Click(object sender, RoutedEventArgs e)
        {
            calculatorManager.SelectedSystem = new CosSystemEquations();
            TextBoxL.IsEnabled = false;
            TextBoxL.Text = "0,5";
        }

        private void ButtonSystemExp_Click(object sender, RoutedEventArgs e)
        {
            calculatorManager.SelectedSystem = new ExpSystemEquations();
            TextBoxL.IsEnabled = true;
            TextBoxL.Text = "";
        }

        private void MethodComboBox_SelectionChanged(object sender, SelectionChangedEventArgs e)
        {
            if (MethodComboBox.SelectedItem is ComboBoxItem selectedItem)
            {
                string selectedText = selectedItem.Content.ToString();
                if (selectedText == "Newton's Method")
                {
                    calculatorManager.SelectedMethod = new NewtonMethod();
                }
                else if (selectedText == "Secant Method")
                {
                    calculatorManager.SelectedMethod = new SecantMethod();
                }
            }
        }

        private void SolveButton_Click(object sender, RoutedEventArgs e)
        {
            if (double.TryParse(TextBoxA.Text, out double a) &&
                double.TryParse(TextBoxB.Text, out double b) &&
                double.TryParse(TextBoxC.Text, out double c) &&
                double.TryParse(TextBoxD.Text, out double d) &&
                double.TryParse(TextBoxK.Text, out double k) &&
                double.TryParse(TextBoxL.Text, out double l) &&
                double.TryParse(TextBoxM.Text, out double m) &&
                double.TryParse(TextBoxN.Text, out double n) &&
                double.TryParse(TextBoxQ.Text, out double q) &&
                double.TryParse(TextBoxR.Text, out double r) &&
                double.TryParse(TextBoxX0.Text, out double x0) &&
                double.TryParse(TextBoxY0.Text, out double y0) &&
                int.TryParse(TextBoxEccuracy.Text, out int e_start))
            {
                double E = 1.0 / Math.Pow(10, e_start);

                calculatorManager.InitializeSystem(a, b, c, d, k, l, m, n, q, r);
                calculatorManager.SetInitialValues(x0, y0, E);

                if (calculatorManager.Solve())
                {
                    labelX.Content = string.Format($"x = {{0:F{e_start}}}", calculatorManager.X);
                    labelY.Content = string.Format($"y = {{0:F{e_start}}}", calculatorManager.Y);
                    labelPC.Content = $"Practical complexity: {calculatorManager.IterationCounter}";
                }
                else
                {
                    labelX.Content = "x = Ø";
                    labelY.Content = "y = Ø";
                    labelPC.Content = $"Practical complexity: {calculatorManager.IterationCounter}";
                }
            }
            else
            {
                MessageBox.Show("Please enter valid numbers in all Fields of Coefficients.");
            }
        }

        private void DownloadFileButton_Click(object sender, RoutedEventArgs e)
        {
            if (calculatorManager.SelectedSystem == null || calculatorManager.SelectedMethod == null)
            {
                MessageBox.Show("Please select system and method first.");
                return;
            }

            FileWorker fileHandler = new FileWorker(calculatorManager);
            string fileContent = fileHandler.GenerateFileContent();
            fileHandler.SaveFile(fileContent);
        }
    }

    public interface ISystemEquations
    {
        void Initialize(double a, double b, double c, double d, double k, double l, double m, double n, double q, double r);
        double Equation1(double x, double y);
        double Equation2(double x, double y);
    }

    public class PowSystemEquations : ISystemEquations
    {
        internal double A, B, C, D, K, L, M, N, Q, R;

        public void Initialize(double a, double b, double c, double d, double k, double l, double m, double n, double q, double r)
        {
            A = a;
            B = b;
            C = c;
            D = d;
            K = k;
            L = l;
            M = m;
            N = n;
            Q = q;
            R = r;
        }

        public double Equation1(double x, double y)
        {
            return A * Math.Pow(x, K) + B * Math.Pow(y, N) - R;
        }

        public double Equation2(double x, double y)
        {
            return C * Math.Pow(x, L) - D * Math.Pow(y, M) - Q;
        }
    }

    public class CosSystemEquations : ISystemEquations
    {
        internal double A, B, C, D, K, L, M, N, Q, R;

        public void Initialize(double a, double b, double c, double d, double k, double l, double m, double n, double q, double r)
        {
            A = a;
            B = b;
            C = c;
            D = d;
            K = k;
            L = l;
            M = m;
            N = n;
            Q = q;
            R = r;
        }

        public double Equation1(double x, double y)
        {
            return A * Math.Pow(Math.Cos(x), K) + B * y - R;
        }

        public double Equation2(double x, double y)
        {
            return C * Math.Sqrt(x) - D * Math.Pow(y, M) - Q;
        }
    }

    public class ExpSystemEquations : ISystemEquations
    {
        internal double A, B, C, D, K, L, M, N, Q, R;

        public void Initialize(double a, double b, double c, double d, double k, double l, double m, double n, double q, double r)
        {
            A = a;
            B = b;
            C = c;
            D = d;
            K = k;
            L = l;
            M = m;
            N = n;
            Q = q;
            R = r;
        }

        public double Equation1(double x, double y)
        {
            return Math.Exp(A * Math.Pow(x, K)) + B * Math.Pow(y, N) - R;
        }

        public double Equation2(double x, double y)
        {
            return Math.Exp(C * Math.Pow(x, L)) - D * Math.Pow(y, M) - Q;
        }
    }

    public interface IMethod
    {
        void Solve(ISystemEquations system, double x0, double y0, double e, out double x, out double y, out int iterations);
    }

    public class NewtonMethod : IMethod
    {
        public void Solve(ISystemEquations system, double x0, double y0, double e, out double x, out double y, out int iterations)
        {
            double tol = 1e-6;
            int maxIter = 10000;
            x = x0;
            y = y0;

            iterations = 0;

            for (int i = 0; i < maxIter; i++)
            {
                iterations++;

                double f1Val = system.Equation1(x, y);
                double f2Val = system.Equation2(x, y);

                if (Math.Abs(f1Val) < tol && Math.Abs(f2Val) < tol)
                {
                    return;
                }

                double J11 = (system.Equation1(x + tol, y) - f1Val) / tol;
                double J12 = (system.Equation1(x, y + tol) - f1Val) / tol;
                double J21 = (system.Equation2(x + tol, y) - f2Val) / tol;
                double J22 = (system.Equation2(x, y + tol) - f2Val) / tol;

                double detJ = J11 * J22 - J12 * J21;

                if (Math.Abs(detJ) < tol)
                {
                    x = double.NaN;
                    y = double.NaN;
                    MessageBox.Show("The determinant of the Jacobian is very small. Root search stopped.");
                    return;
                }

                double dx = (f1Val * J22 - f2Val * J12) / detJ;
                double dy = (f2Val * J11 - f1Val * J21) / detJ;

                double newX = x - dx;
                double newY = y - dy;

                if (Math.Abs(newX - x) < e && Math.Abs(newY - y) < e)
                {
                    x = newX;
                    y = newY;
                    return;
                }

                x = newX;
                y = newY;
            }
            x = double.NaN;
            y = double.NaN;
            MessageBox.Show("Newton's method did not converge.");
        }
    }

    public class SecantMethod : IMethod
    {
        public void Solve(ISystemEquations system, double x0, double y0, double e, out double x, out double y, out int iterations)
        {
            double tol = 1e-6;
            int maxIter = 10000;
            x = x0;
            y = y0;

            iterations = 0;

            double xPrev = x + 0.1;
            double yPrev = y + 0.1;

            double[,] J = new double[2, 2];
            J[0, 0] = (system.Equation1(x, y) - system.Equation1(xPrev, y)) / (x - xPrev);
            J[0, 1] = (system.Equation1(x, y) - system.Equation1(x, yPrev)) / (y - yPrev);
            J[1, 0] = (system.Equation2(x, y) - system.Equation2(xPrev, y)) / (x - xPrev);
            J[1, 1] = (system.Equation2(x, y) - system.Equation2(x, yPrev)) / (y - yPrev);

            double[] F = { system.Equation1(x, y), system.Equation2(x, y) };

            for (int i = 0; i < maxIter; i++)
            {
                iterations++;

                double f1Val = system.Equation1(x, y);
                double f2Val = system.Equation2(x, y);

                if (Math.Abs(f1Val) < tol && Math.Abs(f2Val) < tol)
                {
                    return;
                }

                double detJ = J[0, 0] * J[1, 1] - J[0, 1] * J[1, 0];

                if (Math.Abs(detJ) < tol)
                {
                    x = double.NaN;
                    y = double.NaN;
                    MessageBox.Show("The determinant of the Jacobian is very small. Root search stopped.");
                    return;
                }

                double dx = (f1Val * J[1, 1] - f2Val * J[0, 1]) / detJ;
                double dy = (f2Val * J[0, 0] - f1Val * J[1, 0]) / detJ;

                double newX = x - dx;
                double newY = y - dy;

                if (Math.Abs(newX - x) < e && Math.Abs(newY - y) < e)
                {
                    x = newX;
                    y = newY;
                    return;
                }

                double[] deltaX = { newX - x, newY - y };
                double[] newF = { system.Equation1(newX, newY), system.Equation2(newX, newY) };
                double[] deltaF = { newF[0] - F[0], newF[1] - F[1] };

                double[] JDeltaX = { J[0, 0] * deltaX[0] + J[0, 1] * deltaX[1], J[1, 0] * deltaX[0] + J[1, 1] * deltaX[1] };
                double factor = 1.0 / (deltaX[0] * deltaX[0] + deltaX[1] * deltaX[1]);
                double[] correction = { deltaF[0] - JDeltaX[0], deltaF[1] - JDeltaX[1] };

                for (int k = 0; k < 2; k++)
                {
                    for (int j = 0; j < 2; j++)
                    {
                        J[k, j] += factor * correction[k] * deltaX[j];
                    }
                }

                xPrev = x;
                yPrev = y;
                x = newX;
                y = newY;
                F[0] = newF[0];
                F[1] = newF[1];
            }
            x = double.NaN;
            y = double.NaN;
            MessageBox.Show("Secant's method did not converge.");
        }
    }

    public class CalculatorManager
    {
        public ISystemEquations SelectedSystem { get; set; }
        public IMethod SelectedMethod { get; set; }
        public double X { get; private set; }
        public double Y { get; private set; }
        public double X0 { get; private set; }
        public double Y0 { get; private set; }
        public double E { get; private set; }
        public int IterationCounter { get; private set; }

        public void InitializeSystem(double a, double b, double c, double d, double k, double l, double m, double n, double q, double r)
        {
            SelectedSystem.Initialize(a, b, c, d, k, l, m, n, q, r);
        }

        public void SetInitialValues(double x0, double y0, double e)
        {
            X0 = x0;
            Y0 = y0;
            E = e;
        }

        public bool Solve()
        {
            if (SelectedSystem == null || SelectedMethod == null)
            {
                MessageBox.Show("System and Method must be selected before solving.");
            }

            double resultX, resultY;
            int resultIterationCounter;
            SelectedMethod.Solve(SelectedSystem, X0, Y0, E, out resultX, out resultY, out resultIterationCounter);
            X = resultX;
            Y = resultY;
            IterationCounter = resultIterationCounter;

            return !double.IsNaN(X) && !double.IsNaN(Y);
        }
    }

    public class FileWorker
    {
        private CalculatorManager calculatorManager;

        public FileWorker(CalculatorManager calculatorManager)
        {
            this.calculatorManager = calculatorManager;
        }

        public string GenerateFileContent()
        {
            string systemCondition = "";
            if (calculatorManager.SelectedSystem is PowSystemEquations)
            {
                var system = calculatorManager.SelectedSystem as PowSystemEquations;
                systemCondition = $"System: \n" +
                                  $"Equation 1: Ax^{system.K} + By^{system.N} = R\n" +
                                  $"Equation 2: Cx^{system.L} - Dy^{system.M} = Q\n\n";
            }
            else if (calculatorManager.SelectedSystem is CosSystemEquations)
            {
                var system = calculatorManager.SelectedSystem as CosSystemEquations;
                systemCondition = $"System: \n" +
                                  $"Equation 1: {system.A} * Cos(x)^{system.K} + {system.B} * y = {system.R}\n" +
                                  $"Equation 2: {system.C} * Sqrt(x) - {system.D} * y^{system.M} = {system.Q}\n\n";
            }
            else if (calculatorManager.SelectedSystem is ExpSystemEquations)
            {
                var system = calculatorManager.SelectedSystem as ExpSystemEquations;
                systemCondition = $"System: \n" +
                                  $"Equation 1: Exp({system.A} * x^{system.K}) + {system.B} * y^{system.N} = {system.R}\n" +
                                  $"Equation 2: Exp({system.C} * x^{system.L}) - {system.D} * y^{system.M} = {system.Q}\n\n";
            }

            string methodCondition = "";
            if (calculatorManager.SelectedMethod is NewtonMethod)
            {
                methodCondition = $"Method: Newton's Method\n\n";
            }
            else if (calculatorManager.SelectedMethod is SecantMethod)
            {
                methodCondition = $"Method: Secant Method\n\n";
            }


            string initialValues = $"Initial approximations:\n" +
                                   $"x0 = {calculatorManager.X0}\n" +
                                   $"y0 = {calculatorManager.Y0}\n\n";

            string accuracy = $"Accuracy: e = {-Math.Log10(calculatorManager.E)}\n\n";

            string formattedX = string.Format($"{{0:F{-Math.Log10(calculatorManager.E)}}}", calculatorManager.X);
            string formattedY = string.Format($"{{0:F{-Math.Log10(calculatorManager.E)}}}", calculatorManager.Y);

            string solution = $"Solution:\n" +
                              $"x = {formattedX}\n" +
                              $"y = {formattedY}\n" +
                              $"Practical complexity: {calculatorManager.IterationCounter}";

            return systemCondition + methodCondition + initialValues + accuracy + solution;
        }

        public void SaveFile(string content)
        {
            SaveFileDialog saveFileDialog = new SaveFileDialog
            {
                Filter = "Text file (*.txt)|*.txt",
                Title = "Save Solution",
                FileName = "Solution.txt"
            };

            if (saveFileDialog.ShowDialog() == true)
            {
                File.WriteAllText(saveFileDialog.FileName, content);
            }
        }
    }
}



