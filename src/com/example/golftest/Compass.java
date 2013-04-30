package com.example.golftest;

import java.util.Random;

import android.app.Activity; 
import android.os.Bundle; 
import android.os.Handler;
import android.view.Gravity;
import android.widget.LinearLayout;
import android.widget.RelativeLayout;
import android.widget.TextView;

public class Compass extends Activity{ 

	private TextView t;
	private ImageCompass imageCompass;
	private AudioService audioService;
	private Random random = new Random();

	private final Runnable mUpdateUITimerTask = new Runnable() {
		public void run() {
			while(true){
				receiveMyMessage();
				try {
					Thread.sleep(1000);
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
	};

	public static Thread performOnBackgroundThread(final Runnable runnable) {
		final Thread t = new Thread() {
			@Override
			public void run() {
				try {
					runnable.run();
				} finally {

				}
			}
		};
		t.start();
		return t;
	}

	private Handler mHandler = new Handler();

	public void receiveMyMessage() {
		final float val = random.nextFloat();
		final float degree = val * 180;

		mHandler.post(new Runnable() {
			@Override
			public void run() {
				if( audioService != null ) {
					audioService.doUpdate( System.currentTimeMillis() );
					// This gets executed on the UI thread so it can safely modify Views
					t.setText( Float.toString( audioService.peakSpecRead() ) );
					/* set the degree of the compass */
					imageCompass.setDirection( (int)degree );
				}
			}
		});
	}
	
	@Override
	protected void onCreate(Bundle icicle) { 
		super.onCreate(icicle); 

		/*
		 * manually define the layout
		 */
		LinearLayout layout = new LinearLayout(this);
        layout.setGravity(Gravity.CENTER);
        layout.setOrientation(LinearLayout.VERTICAL);

        /*
         * manually create the text view and position it
         */
		t = new TextView(this);
		RelativeLayout.LayoutParams params = new RelativeLayout.LayoutParams(100, 80);
		params.leftMargin = 50;
		params.topMargin = 60;
		t.setWidth(200);
		t.setHeight(100);
	    t.setText("Step One: blast egg");
	    layout.addView(t, params);

	    /*
	     * manually create the imageview with compass picture and position it
	     */
	    RelativeLayout.LayoutParams imageViewParams = new RelativeLayout.LayoutParams(100, 100);
		params.leftMargin = 150;
		params.topMargin = 160;
		imageCompass = new ImageCompass(this);
	    layout.addView(imageCompass, imageViewParams);

	    /*
	     * set the content view to this layout
	     */
	    setContentView(layout);
	    
	    /**
	     * set up the audio service
	     */
        audioService = new AudioService();
	    audioService.measureStart();
	    /*
	     * wake up the thread to keep updating values every second
	     */
		performOnBackgroundThread(mUpdateUITimerTask);
	}
	
    @Override
    protected void onStop() {
        super.onStop();
        if(audioService!=null) {
        	audioService.onDestroy();
        }
    }	
}