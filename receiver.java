package TermProject;
import java.io.*;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;

public class receiver
{
	static int aI = 0; // global variable used to keep track of the index of args
	static String[] input;
	static double[][] data, J;
	static double[] X0, B, XV, Sk;
	static double TV;
	static int M;
	static double c = Double.parseDouble("2.997924580000000000E+08");
	static double R = Double.parseDouble("6.367444500000000000E+06");
	static double S = Double.parseDouble("8.616408999999999651E+04");
	static ArrayList<String> vals;
	static ArrayList<String> output;
	static NumberFormat formatter;
	
	public static void main(String[] args) throws IOException
	{
		formatter = new DecimalFormat("#0.00");
		
		
		output = new ArrayList<String>();
		vals = new ArrayList<String>();
		ReadInput();
		X0 = new double[]{ 0, 0, 0 };
		while(aI<input.length)
		{
			data = ReadData(); // data is an m by n matrix
			ComputeLocation();
		}
		WriteOutput();
	}
	
	static void WriteOutput()
	{
		try(Writer writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream("receiver.log"),"utf-8")))
		{
			for(int i = 0; i<output.size(); i++) { writer.write(output.get(i)+"\n"); }
			writer.close();
		}
		catch(IOException e) { System.out.println("Something went wrong."); }
	}
	
	static void ReadInput() throws IOException
	{
		BufferedReader In = new BufferedReader(new InputStreamReader(System.in));
		String lines;
		while((lines = In.readLine()) != null)
		{
			for(String st : lines.split(" ")) { vals.add(st); }
		}
		int cap = vals.size();
		input = new String[cap];
		for(int i = 0; i < cap; i++) { input[i] = vals.get(i); }
		In.close();
	}
	
	
	static void ComputeLocation()
	{
		XV = NewtonsMethod(X0, 0);
		TV = TimeAt(XV);
		String[] retval = PositionInGeographic(new String[]{""+TV, ""+XV[0], ""+XV[1], ""+XV[2]});
		String ret = "";
		for(int i = 0; i<10; i++)
		{
			ret = ret + formatter.format(Double.parseDouble(retval[i])) + " ";
		}
		System.out.println(ret);
		output.add(ret);
	}
	
	static double[] NewtonsMethod(double[] x, int d)
	{
		Sk = Solve3x3(Jacobian(x), gradf(x));
		double[] retval = x.clone();
		for(int i = 0; i<3; i++) { retval[i] -= Sk[i]; }
		if(Math.abs(diff(retval,x)[0])<0.0000001d && Math.abs(diff(retval,x)[1])<0.0000001d && Math.abs(diff(retval,x)[2])<0.0000001d) { return retval; }
		else if(d>9) { return null; }
		else
		{
			d++;
			return NewtonsMethod(retval, d);
		}
	}
	
	static double[][] Jacobian(double[] x)
	{
		double[][] retval = new double[3][3];
		for(int i = 0; i < 3; i++)
		{
			for(int j = 0; j < 3; j++)
			{
				retval[i][j] = J_ij(i, j, x);
			}
		}
		return retval;
	}
	
	static double twonorm(double[] a)
	{
		double sum = 0.0d;
		for(Double d : a){ sum = sum + d*d; }
		return Math.sqrt(sum);
	}
	
	static double[] diff(double[] a, double[] b)
	{
		if(a.length==b.length)
		{
			double[] retval = new double[a.length];
			for(int i = 0; i<a.length; i++) { retval[i]=a[i]-b[i]; }
			return retval;
		}
		return null;
	}
	
	static double[][] ReadData()
	{
		int bounds = CalculateBounds();
		double[][] retval = new double[bounds][4];
		// while the current first input is greater than the previous first input
		for(int j = 0; j < bounds; j++)
		{
			// Read 5 values, input[aI] through input[aI+4]
			// Input 1 is vehicle id - Don't even need this? Just throw it away
			// Input 2 is the time at which the vehicle broadcasts
			retval[j][0] = Double.parseDouble(input[aI + 1]);
			// Input 3, 4, 5 are the Cartesian coordinates (x,y,z) of the satellite at time of broadcast
			retval[j][1] = Double.parseDouble(input[aI + 2]);
			retval[j][2] = Double.parseDouble(input[aI + 3]);
			retval[j][3] = Double.parseDouble(input[aI + 4]);


			// Increment aI by 5
			aI = aI + 5;
		}
		M = bounds;
		return retval;
	}
	
	// Calculate the bounds of each chunk of satellite data.
	static int CalculateBounds()
	{
		int m = 1;
		while(aI+1+5*m < input.length && Math.abs(Double.parseDouble(input[aI+1+5*(m)])-Double.parseDouble(input[aI+1+5*(m-1)]))<0.5d) { m++; }
		return m;
	}
	
	// 3x3 Linear System solver. (Will return NaN NaN NaN if A is singular or there are other elimination problems)
	static double[] Solve3x3(double[][] A, double[] b)
	{
		double x = (A[0][1]*A[1][2]*b[2] - A[0][1]*A[2][2]*b[1] - A[0][2]*A[1][1]*b[2] + A[0][2]*A[2][1]*b[1] + A[1][1]*A[2][2]*b[0] - A[1][2]*A[2][1]*b[0])/(A[0][0]*A[1][1]*A[2][2] - A[0][0]*A[1][2]*A[2][1] - A[0][1]*A[1][0]*A[2][2] + A[0][1]*A[1][2]*A[2][0] + A[0][2]*A[1][0]*A[2][1] - A[0][2]*A[1][1]*A[2][0]);
		double y = -1 * (A[0][0]*A[1][2]*b[2] - A[0][0]*A[2][2]*b[1] - A[0][2]*A[1][0]*b[2] + A[0][2]*A[2][0]*b[1] + A[1][0]*A[2][2]*b[0] - A[1][2]*A[2][0]*b[0])/(A[0][0]*A[1][1]*A[2][2] - A[0][0]*A[1][2]*A[2][1] - A[0][1]*A[1][0]*A[2][2] + A[0][1]*A[1][2]*A[2][0] + A[0][2]*A[1][0]*A[2][1] - A[0][2]*A[1][1]*A[2][0]);
		double z = (A[0][0]*A[1][1]*b[2] - A[0][0]*A[2][1]*b[1] - A[0][1]*A[1][0]*b[2] + A[0][1]*A[2][0]*b[1] + A[1][0]*A[2][1]*b[0] - A[1][1]*A[2][0]*b[0])/(A[0][0]*A[1][1]*A[2][2] - A[0][0]*A[1][2]*A[2][1] - A[0][1]*A[1][0]*A[2][2] + A[0][1]*A[1][2]*A[2][0] + A[0][2]*A[1][0]*A[2][1] - A[0][2]*A[1][1]*A[2][0]);
		return new double[]{x,y,z};
	}
	
	// Returns the (i, j)-th component of the square symmetric Jacobian for the least squares problem
	static double J_ij(int i, int j, double[] x)
	{
		double sum = 0.0d;
		for(int k = 0; k<M-1; k++) { sum += gradF_ij(k, i, x) * gradF_ij(k, j, x); }
		return 2.0d * sum;
	}
	
	// Returns the (i, j)-th component of the Gradient of F(x)
	static double gradF_ij(int i, int j, double[] x)
	{
		double retval = 0.5d * (2.0d * x[j] - 2.0d * data[i][j+1]) / (Math.sqrt(Math.pow((x[0]-data[i][1]),2) + Math.pow((x[1]-data[i][2]),2) + Math.pow((x[2]-data[i][3]),2)));
		retval -= 0.5d * (2.0d * x[j] - 2.0d * data[i+1][j+1]) / (Math.sqrt(Math.pow((x[0]-data[i+1][1]),2) + Math.pow((x[1]-data[i+1][2]),2) + Math.pow((x[2]-data[i+1][3]),2)));
		return retval;
	}
	
	static double[] gradf(double[] x)
	{
		double[][] diffs = new double[M][3];
		for(int j = 0; j<M; j++) { diffs[j] = diff(new double[]{data[j][1],data[j][2],data[j][3]},x); }
		
		double[] N = new double[M];
		for(int j = 0; j<M; j++) { N[j] = twonorm(diffs[j]); }
		
		double[] A = new double[M-1];
		for(int j = 0; j<M-1; j++) { A[j] = N[j+1]-N[j]-c*(data[j][0]-data[j+1][0]); }
		
		double[][] XYZ = new double[3][M-1];
		for(int i = 0; i<3; i++)
		{
			for(int j = 0; j<M-1; j++)
			{
				XYZ[i][j] = diffs[j][i]/N[j] - diffs[j+1][i]/N[j+1];
			}
		}
		
		double[] retval = new double[3];
		for(int j = 0; j<3; j++)
		{
			retval[j] = 0.0d;
			for(int i = 0; i<M-1; i++) { retval[j] += A[i]*XYZ[j][i]; }
			retval[j] = 2.0d*retval[j];
		}
		return retval;
	}
	
	// Once Xv is computed, this returns the corresponding time Tv
	static double TimeAt(double[] Xv)
	{
		double retval = data[0][0] + (1/c) * Math.sqrt(Math.pow(Xv[0]-data[0][1], 2)+Math.pow(Xv[1]-data[0][2], 2)+Math.pow(Xv[2]-data[0][3], 2));
		return retval;
	}
	
	private static String[] PositionInGeographic(String[] a)
	{
		double t = Double.parseDouble(a[0]);
		double x1 = Double.parseDouble(a[1]);
		double y1 = Double.parseDouble(a[2]);
		double z1 = Double.parseDouble(a[3]);
		
		double[] xyz = R3(-2*Math.PI*t/S, new double[]{x1,y1,z1});
		double x = xyz[0];
		double y = xyz[1];
		double z = xyz[2];
		
		double psi;
		if(x*x+y*y == 0) {
			if(z>=0) { psi = Math.PI/2;	}
			else { psi = -1*Math.PI/2;	}
		}
		else { psi = Math.atan2(z,Math.sqrt(x*x+y*y)); 
		}
		double lambda;
		if(x>0 && y>0) { lambda = Math.atan2(y,x); }
		else if(x < 0) { lambda = Math.PI + Math.atan2(y,x); }
		else { lambda = 2*Math.PI + Math.atan2(y,x); }
		lambda-=Math.PI;
		
		String[] Psi = RadiansToDegrees(psi);
		String[] Lambda = RadiansToDegrees(lambda);
		
		double h = Math.sqrt(x*x + y*y + z*z) - R;
		
		return new String[]{""+t, Psi[0], Psi[1], Psi[2], Psi[3], Lambda[0], Lambda[1], Lambda[2], Lambda[3], ""+h};
	}
	private static String[] RadiansToDegrees(double a)
	{
		double 	b  		= 	a*180/Math.PI;
		if(a<0){b 		= 	-1*b;}
		int 	d 		= 	((int) Math.floor(b));
		int 	m		= 	((int) Math.floor(60*(b-d)));
		double 	s 		= 	60*(60*(b-d)-m);
		if(a<0)	{ return new String[]{""+d,""+m,""+s, ""+(-1)};	}
		return new String[]{""+d,""+m,""+s, ""+1};
	}
	private static double[] R3(double alpha, double[] x)
	{
		return new double[] { Math.cos(alpha)*x[0] - Math.sin(alpha)*x[1], Math.sin(alpha)*x[0] + Math.cos(alpha)*x[1], x[2] };
	}
}
