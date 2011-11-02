package circdesigna.TripleSim;

import java.util.LinkedList;
import java.util.Queue;

import circdesigna.TripleSim.ReactionGraph3X.GraphNode;

public class PulseEvents {
	public static class Event implements Comparable<Event>{
		public double t;
		public double amount;

		public Event(double t, double amount){
			this.t = t;
			this.amount = amount;
		}

		public int compareTo(Event o) {
			if (t < o.t){
				return -1;
			}
			if (t > o.t){
				return 1;
			}
			return 0;
		}
	}
	public Event[] schedule;
	private Queue<Event> events = new LinkedList<Event>();
	private GraphNode v;
	public PulseEvents(GraphNode init, double ... events) {
		v = init;
		schedule = new Event[events.length/2];
		for(int i = 0; i < events.length/2; i++){
			schedule[i] = new Event(events[i*2],events[i*2+1]);
		}
	}
	public void reset(){
		events.clear();
		for(Event q : schedule){
			events.offer(q);
		}
	}

	public double[] handlePulses(double[] concs, double time, double step){
		boolean shocked = false;
		while (!events.isEmpty() && time >= events.peek().t){
			shocked = true;
			Event poll = events.poll();
			time = poll.t;
			concs[v.index] += poll.amount;
		}
		return new double[]{time, shocked?Math.min(step,1e-4):step};
	}
}
