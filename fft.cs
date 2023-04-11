using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Drawing;

namespace DiaGui_2
{
    /// <summary>
    /// FFT code based on Cooley-Tukey algorithm
    /// 
    /// </summary>
    class FFT
    {
        private int M, N, Fs, currentIndex;
      	private double[] cos;
        private double[] sin;

        private double[] real;
        private double[] imag;
        private double[] magnitude;
        private int[] intMag;
        private double[] avgmag;

        private float[] sampleData;
        private int dataCounter = 0;
        private float maxFreq = 0;
        //private static bool dataCollected = false;
        
        public FFT(int N, int Fs)
        {
            this.N = N;
            this.Fs = Fs;
            M = (int) (Math.Log(N) / Math.Log(2));
            //--- Make sure N is a power of 2
		    if (N != (1 << M)) {
                throw new Exception("FFT length must be a power of 2");
            }
            real = new double[N];
			imag = new double[N];
			magnitude = new double[N];
            intMag = new int[N];
			avgmag = new double[N];
            sampleData = new float[N];

	      	//initialize avgmag array
			for (int xx=0; xx<N; xx++){
				avgmag[xx] = 0;
			}
	      
	      	//pre-compute tables
			cos = new double[N / 2];
			sin = new double[N / 2];

			for (int i = 0; i < N / 2; i++) {
				cos[i] = Math.Cos(-2 * Math.PI * i / N);
				sin[i] = Math.Sin(-2 * Math.PI * i / N);
    		}

	    }
        public void PutSample(float newSample)
        {
            int n2 = N/2;
            //--- update sample
            sampleData[currentIndex] = newSample;
            //real[currentIndex] = (double) newSample;
            ShiftSamples();

            //--- update magnitude
            magnitude[currentIndex] = GetMagnitude(real[currentIndex], imag[currentIndex]);

            //--- get the frequency with largest amplitude
            maxFreq = GetFrequency(getMaxFreqIndex());
            
            currentIndex++;
            if (currentIndex == N){// / 256) {
                PutAll(sampleData);
                //--- get FFT
                System.Threading.Thread fftThread = new System.Threading.Thread(GetFft);
                fftThread.Priority = System.Threading.ThreadPriority.AboveNormal;
                fftThread.Start();
                
            }

            if (currentIndex == N)
            {
                //PutAll(sampleData);
                currentIndex = 0;
            }
            dataCounter++;
        }
        public void PutSingleSample(float newSample)
        {
            //--- shift
            Array.Copy(real, 1, real, 0, N - 1);
            Array.Copy(imag, 1, imag, 0, N - 1);
            
            //--- insert new sample
            real[N - 1] = (double)newSample;
            imag[N - 1] = 0;

            //--- update magnitude
            magnitude[currentIndex] = GetMagnitude(real[currentIndex], imag[currentIndex]);

            //--- FFT
            dataCounter++;
            if (dataCounter == N / 64)
            {
                dataCounter = 0;
                //--- get FFT
                System.Threading.Thread fftThread = new System.Threading.Thread(GetFft);
                fftThread.Priority = System.Threading.ThreadPriority.AboveNormal;
                fftThread.Start();
                
            }
            currentIndex++;
            if (currentIndex == N)
            {
                currentIndex = 0;
            }



        }
		public void PutAll(float[] inData){
			for (int i=0; i<N; i++){
				real[i] = (double) inData[i];
				imag[i] = 0;
			}
		}
        public void GetFft()
        {

            int i, j, k, n1, n2, a;
        	double c, s, t1, t2;

	        // clear imag
	        for (i=0; i < N; i++){
		        imag[i] = 0;
	        }
	        // Bit-reverse
	        j = 0;
	        n2 = N / 2;
	        for (i = 1; i < N - 1; i++) {
		        n1 = n2;
		        while (j >= n1) {
			        j = j - n1;
			        n1 = n1 / 2;
		        }
		        j = j + n1;

		        if (i < j) {
                    t1 = real[i];
                    real[i] = real[j];
                    real[j] = t1;
			        t1 = imag[i];
                    imag[i] = imag[j];
                    imag[j] = t1;
		        }
	        }

	        // FFT
	        n1 = 0;
	        n2 = 1;

	        for (i = 0; i < M; i++) {
		        n1 = n2;
		        n2 = n2 + n2;
		        a = 0;

		        for (j = 0; j < n1; j++) {
			        c = cos[a];
			        s = sin[a];
			        a += 1 << (M - i - 1);

			        for (k = j; k < N; k = k + n2) {
				        t1 = c * real[k + n1] - s * imag[k + n1];
				        t2 = s * real[k + n1] + c * imag[k + n1];
				        real[k + n1] = real[k] - t1;
				        imag[k + n1] = imag[k] - t2;
				        real[k] = real[k] + t1;
				        imag[k] = imag[k] + t2;
			        }
		        }
	        }
        }
        public double GetMagnitude(double real, double imag)
        {
	        return 4 * (Math.Sqrt(real * real + imag * imag))/N;
        }
        public void GetMagnitudeAll(int endFreq)
        {
	        int i;
	        int endIndex;
	        endIndex = endFreq * N / Fs;
	        for (i=0; i<endIndex; i++){
                magnitude[i] = GetMagnitude(real[i], imag[i]);
                intMag[i] = (int)magnitude[i];
	        }
        }
        public void Draw(Graphics g, int top, double scale, int endFreq, bool hasCursor, Color color)
        {
            Pen pen;
            //pen = new System.Drawing.Pen(System.Drawing.Color.Orange);
            pen = new System.Drawing.Pen(color);
            pen.Width = 2;
            int freq = 20;
            SolidBrush brush = new SolidBrush(Color.Yellow);
            Font font = new Font("Courier", 6);
            g.SmoothingMode = System.Drawing.Drawing2D.SmoothingMode.AntiAlias;

            //--- draw max frequency value
            g.DrawString(maxFreq.ToString("0.0 Hz"), font, brush, 200, 30);

            int endIndex;
            endIndex = endFreq * N / Fs;

            //--- draw
            int i = 0;
            int j = 0;
            int k = 0;
            for (i = 0; i < endIndex; i++)
            {
                g.DrawLine(pen, i, top, i, top - (int) (magnitude[i] * scale));

                //--- draw x-axis grid lines
                j++;
                k++;
                if (j == 21)        // 10 Hz interval
                {
                    pen.Color = Color.Yellow;
                    
                    g.DrawLine(pen, i, top, i, top - 5);
                    j = 0;
                    //--- draw baseline
                    g.DrawLine(pen, 0, top, i, top);
                    //--- draw freq values
                    g.DrawString(freq.ToString(), font, brush, i - 5, top + 2);

                    freq += 20;
                    if (freq > 240)
                    {
                        freq = 20;
                    }

                    pen.Color = color;
                }

            }
            /*
            if (isActive)
            {
                if (dataCounter < N)
                {
                    Font font = new Font("Courier", 8);
                    SolidBrush brush = new SolidBrush(Color.Orange);
                    g.DrawString("collecting data...", font, brush, 80, top - 50);
                }
                else
                {
                    dataCollected = true;
                }
                Console.WriteLine("" + dataCounter);
            }
            */
        }
        public float GetFrequency(int index)
        {
            float freq = 1.0f * index * Fs / N;
            return freq;
        }

        public bool isActive { get; set; }
        private int getMaxFreqIndex()
        {
            double max = 0;
            //double maxFreq = magnitude.Max();
            for (int i = 0; i < N / 2; i++)
            {
                if (magnitude[i] > max)
                {
                    max = magnitude[i];
                }
            }

            return Array.IndexOf(magnitude, max);
        }
        public void setMaxFreq(int maxFreq)
        {
            this.maxFreq = maxFreq;
        }
        private void ShiftSamples()
        {

        }
    }
}
