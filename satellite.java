package TermProject;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Scanner;

public class satellite
{
	static double Pi, c, R, S;
	static double[][] V;
	static double[] Xv, Xs;
	static double Tv, Ts;
	static String path;
	static ArrayList<String> output;
	
	public static void main(String[] args) throws IOException
	{
		path = ".\\";
		
		Ts = 0.0d;
		Xs = null;
		Initialize();
		ReadInput();
		WriteOutput();
	}
	
	private static void ReadInput() throws IOException
	{
		BufferedReader In = new BufferedReader(new InputStreamReader(System.in));
		String s;
		while((s = In.readLine()) != null) { ProcessLine(s); }
	}

	private static void ProcessLine(String s)
	{
		String[] split = s.split(" ");
		double[] xyz = PositionInCartesian(split);
		Tv = Double.parseDouble(split[0]);
		Xv = xyz;
		boolean[] B = HorizonCheck(xyz);
		for(int i = 0; i < B.length; i++)
		{
			if(B[i])
			{
				Ts = ComputeTs(i);
				Xs = Xs(i, Ts);
				String ret = i+" "+Ts+" "+Xs[0]+" "+Xs[1]+" "+Xs[2];
				System.out.println(ret);
				output.add(ret);
				if(Ts>200 && Ts<203)
				{
					System.out.println("Look here!");
					output.add("Look here!");
				}
			}
		}
	}
	
	static void WriteOutput()
	{
		try(Writer writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream("satellite.log"),"utf-8")))
		{
			for(int i = 0; i<output.size(); i++) { writer.write(output.get(i)+"\n"); }
			writer.close();
		}
		catch(IOException e) { System.out.println("Something went wrong."); }
	}
	
	private static double ComputeTs(int v)
	{
		double t0 = Tv;
		double norm = 0.0d;
		for(int i = 0; i < 3; i++)
		{
			norm += Math.pow(Xs(v,Tv)[i],2);
		}
		t0 = t0 - Math.sqrt(norm)/c;
		
		return NewtonsMethod(v,t0,0);
	}
	
	private static void Initialize()throws IOException
	{
		V = new double[24][9];
		output = new ArrayList<String>();
		String temp = "";
		Scanner scan = new Scanner(new File(path+"data.dat")).useDelimiter("\n");
		scan.useDelimiter("\n");
		String[] l = new String[220];
		int counter = 0;
		while(scan.hasNext())
		{
			temp = scan.next();
			l[counter] = temp;
			counter++;
		}
		scan.close();
		Pi = Double.parseDouble(l[0].substring(0, 27));
		c = Double.parseDouble(l[1].substring(0, 27));
		R = Double.parseDouble(l[2].substring(0, 27));
		S = Double.parseDouble(l[3].substring(0, 27));
		
		counter = 4;
		for(int v = 0; v < 24; v++)
		{
			for(int i = 0; i < 9; i++)
			{
				V[v][i] = Double.parseDouble(l[counter].substring(0, 27));
				counter++;
			}
		}
	}
	
	// Computations	
	private static double[] Xs(int v, double t)
	{
		double x = (R+V[v][7])*(V[v][0]*Math.cos((2*Pi*t)/V[v][6]+V[v][8]) + V[v][3]*Math.sin((2*Pi*t)/V[v][6]+V[v][8]));
		double y = (R+V[v][7])*(V[v][1]*Math.cos((2*Pi*t)/V[v][6]+V[v][8]) + V[v][4]*Math.sin((2*Pi*t)/V[v][6]+V[v][8]));
		double z = (R+V[v][7])*(V[v][2]*Math.cos((2*Pi*t)/V[v][6]+V[v][8]) + V[v][5]*Math.sin((2*Pi*t)/V[v][6]+V[v][8]));
		return new double[]{x,y,z,t};
	}
	private static double[] Xs1(int v, double t)
	{
		double x = 2*Pi*(1/V[v][6])*(R+V[v][7])*(V[v][3]*Math.sin((2*Pi*t)/V[v][6]+V[v][8])-V[v][0]*Math.cos((2*Pi*t)/V[v][6]+V[v][8]));
		double y = 2*Pi*(1/V[v][6])*(R+V[v][7])*(V[v][4]*Math.sin((2*Pi*t)/V[v][6]+V[v][8])-V[v][1]*Math.cos((2*Pi*t)/V[v][6]+V[v][8]));
		double z = 2*Pi*(1/V[v][6])*(R+V[v][7])*(V[v][5]*Math.sin((2*Pi*t)/V[v][6]+V[v][8])-V[v][2]*Math.cos((2*Pi*t)/V[v][6]+V[v][8]));
		return new double[]{x,y,z,t};
	}
	private static double f(int v, double t)
	{
		double retval = -1*c*c*(Tv-t)*(Tv-t);
		for(int i = 0; i<3; i++)
		{
			retval += (Xs(v,t)[i]-Xv[i])*(Xs(v,t)[i]-Xv[i]);
		}
		return retval;
	}
	private static double f1(int v, double t)
	{
		double retval = 2*c*c*(Tv-t);
		for(int i = 0; i<3; i++)
		{
			retval += 2*(Xs(v,t)[i]-Xv[i])*Xs1(v,t)[i];
		}
		return retval;
	}
	private static double NewtonsMethod(int v, double Tk, int depth)
	{
		double retval = Tk;
		retval = Tk-(f(v,Tk)/f1(v,Tk));
		if(retval-Tk < 0.01/c) { return retval; }
		else if(depth>=9) { return -1; }
		else { return NewtonsMethod(v, retval, depth+1); }
	}
	
	private static String[] RadiansToDegrees(double a)
	{
		double 	b  		= 	a*180/Pi;
		if(a<0){b 		= 	-1*b;}
		int 	d 		= 	((int) Math.floor(b));
		int 	m		= 	((int) Math.floor(60*(b-d)));
		double 	s 		= 	60*(60*(b-d)-m);
		if(a<0){d		=	-1*d;}
		return new String[]{""+d,""+m,""+s};
	}

	// Utilities
	private static double[] PositionInCartesian(String[] a)
	{
		double 	t 		= 	Double.parseDouble(a[0]);
		double 	theta 	= 	Integer.parseInt(a[4])*DegreesToRadians(new String[]{a[1],a[2],a[3]});
		double 	phi		= 	Integer.parseInt(a[8])*DegreesToRadians(new String[]{a[5],a[6],a[7]});
		double 	h 		= 	Double.parseDouble(a[9]);
		double 	x 		= 	(R + h) * Math.cos(theta) * Math.cos(phi);
		double 	y 		= 	(R + h) * Math.cos(theta) * Math.sin(phi);
		double 	z 		= 	(R + h) * Math.sin(theta);
		double  alpha   =   (2*Pi*t)/S;
		
		return new double[]{Math.cos(alpha)*x-Math.sin(alpha)*y,Math.sin(alpha)*x+Math.cos(alpha)*y,z,t};
	}
	private static double DegreesToRadians(String[] a)
	{
		return 2*Pi*(Integer.parseInt(a[0])/360.0d + Integer.parseInt(a[1])/(360.0d*60) + Double.parseDouble(a[2])/(360*60*60));
	}
	private static boolean[] HorizonCheck(double[] a)
	{
		boolean[] retval = new boolean[24];
		double[] satXYZ;
		for(int i = 0; i < 24; i++)
		{
			satXYZ = Xs(i, a[3]);
			retval[i] = 2*a[0]*(satXYZ[0]-a[0]) + 2*a[1]*(satXYZ[1]-a[1]) + 2*a[2]*(satXYZ[2]-a[2]) > 0;
		}
		return retval;
	}
}
