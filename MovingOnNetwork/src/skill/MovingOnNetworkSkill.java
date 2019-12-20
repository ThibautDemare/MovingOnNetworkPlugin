package skill;

import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.graphstream.algorithm.Dijkstra;
import org.graphstream.graph.Edge;
import org.graphstream.graph.EdgeRejectedException;
import org.graphstream.graph.Graph;
import org.graphstream.graph.IdAlreadyInUseException;
import org.graphstream.graph.Node;
import org.graphstream.graph.Path;
import org.graphstream.graph.implementations.MultiGraph;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.LineString;
import com.vividsolutions.jts.geom.Point;

import msi.gama.common.interfaces.ILocated;
import msi.gama.common.geometry.GeometryUtils;
import msi.gama.metamodel.agent.IAgent;
import msi.gama.metamodel.shape.GamaPoint;
import msi.gama.metamodel.shape.ILocation;
import msi.gama.metamodel.shape.IShape;
import msi.gama.metamodel.topology.ITopology;
import msi.gama.metamodel.topology.filter.In;
import msi.gama.metamodel.topology.graph.GraphTopology;
import msi.gama.precompiler.GamlAnnotations.action;
import msi.gama.precompiler.GamlAnnotations.arg;
import msi.gama.precompiler.GamlAnnotations.doc;
import msi.gama.precompiler.GamlAnnotations.example;
import msi.gama.precompiler.GamlAnnotations.getter;
import msi.gama.precompiler.GamlAnnotations.setter;
import msi.gama.precompiler.GamlAnnotations.skill;
import msi.gama.precompiler.GamlAnnotations.variable;
import msi.gama.precompiler.GamlAnnotations.vars;
import msi.gama.precompiler.GamlAnnotations;
import msi.gama.runtime.GAMA;
import msi.gama.runtime.IScope;
import msi.gama.runtime.exceptions.GamaRuntimeException;
import msi.gama.util.GamaDate;
import msi.gama.util.GamaListFactory;
import msi.gama.util.IList;
import msi.gama.util.graph.GamaGraph;
import msi.gama.util.graph.GraphUtilsGraphStream;
import msi.gama.util.graph.IGraph;
import msi.gama.util.graph._Edge;
import msi.gama.util.graph._Vertex;
import msi.gaml.operators.Cast;
import msi.gaml.skills.Skill;
import msi.gaml.types.IType;

@doc("This skill is intended to move an agent on a network according to speed and length attributes on the edges. When The agent is not already on the graph, we assume that the length is an euclidean length and we use a default speed given by the user.")
@vars({
	@variable (
			name = IKeywordMoNAdditional.DEFAULT_SPEED,
			type = IType.FLOAT,
			init = "19.4444",
			doc = @doc("The speed outside the graph (in meter/second). Default : 70km/h.")),
	@variable(
			name = IKeywordMoNAdditional.PATH_LENGTH,
			type = IType.FLOAT,
			doc = @doc("The length of the computed path.")),
})
@skill(name = IKeywordMoNAdditional.MOVING_ON_NETWORK)
public class MovingOnNetworkSkill extends Skill {

	private class DataSimulation {
		/*
		 * Utils
		 */
		private class DataNetwork {
			private Dijkstra dijkstra;
			private Graph graph;
			private FileSinkDGSFiltered fileSink;
			private GamaGraph gamaGraph;
			private String fileName;
			private String lengthAttribute;
			private String speedAttribute;
		}

		/*
		 * Static variables
		 */
		private int currentCycle = 0;
		private HashMap<String, DataNetwork> listNetwork;
		private DataNetwork currentDataNetwork;
		/*
		 * Non-static variables
		 */
		private double remainingTime;
		private boolean agentFromOutsideToInside;
		private boolean agentInside;
		private boolean agentFromInsideToOutside;
		private boolean agentOnANode;
		private int indexSegment;
		private List<Edge> currentGsPathEdge;
		private List<Node> currentGsPathNode;
		private ILocation currentTarget;
		private double seed;
	}

	// All the data used by this plugin are stored into this variable, itself stored as an attribute of a simulation.
	// Therefore, this plugin can be used with simulation executed in parallel.
	private static DataSimulation currentSimulation;

	/*
	 * Getters and setters
	 */

	public void getCurrentSimulation(final IScope scope) {
		currentSimulation = (DataSimulation) scope.getSimulation().getAttribute("gaml.extensions.moving.on.network.dataSimulation");
		if(currentSimulation == null) {
			currentSimulation = new DataSimulation();
			currentSimulation.seed = scope.getSimulation().getSeed();
			scope.getSimulation().setAttribute("gaml.extensions.moving.on.network.dataSimulation", currentSimulation);
		}
	}

	private ILocation getTarget(final IScope scope) {
		final Object target = scope.getArg(IKeywordMoNAdditional.TARGET, IType.NONE);
		if ( target != null && target instanceof IShape )
			return ((ILocated) target).getLocation();
		return null;
	}

	private ILocation getSource(final IScope scope) {
		final Object source = scope.getArg(IKeywordMoNAdditional.SOURCE, IType.NONE);
		if ( source != null && source instanceof IShape )
			return ((ILocated) source).getLocation();
		return null;
	}

	@getter(IKeywordMoNAdditional.DEFAULT_SPEED)
	public Double getDefaultSpeed(final IAgent agent) {
		return (Double) agent.getAttribute(IKeywordMoNAdditional.DEFAULT_SPEED);
	}

	@setter(IKeywordMoNAdditional.DEFAULT_SPEED)
	public void setDefaultSpeed(final IAgent agent, final double s) {
		agent.setAttribute(IKeywordMoNAdditional.DEFAULT_SPEED, s);
	}

	/*
	 * Actions/methods
	 */

	@action(
			name = "leave_building",
			doc =
			@doc(value = "Inform the building that we leave it.", examples = { @example("do leave_building;") })
			)
	public void leaveBuilding(final IScope scope) throws GamaRuntimeException {
		IAgent v = (IAgent)getCurrentAgent(scope);
		IAgent source = (IAgent) v.getAttribute("source");
		String networkType = (String) v.getAttribute("networkType");
		String destName = source.getName();
		source.setAttribute("lastVehicleDepartureDest_"+networkType+"_"+destName, scope.getClock().getCurrentDate());
	}

	@action(
			name = "add_network",
			args = {
					@arg(name = IKeywordMoNAdditional.NAME, type = IType.STRING, optional = false, doc = @doc("the name of the network.")),
					@arg(name = IKeywordMoNAdditional.GRAPH, type = IType.GRAPH, optional = false, doc = @doc("the network to add.")),
					@arg(name = IKeywordMoNAdditional.LENGTH_ATTRIBUTE, type = IType.STRING, optional = true, doc = @doc("the name of the variable containing the length of an edge.")),
					@arg(name = IKeywordMoNAdditional.SPEED_ATTRIBUTE, type = IType.STRING, optional = true, doc = @doc("the name of the varaible containing the speed on an edge."))
			},
			doc =
			@doc(value = "Add a network.", examples = { @example("do add_network name:'maritime' network:maritime_network length_attribute:'length' speed_attribute:'speed';") })
			)
	public void addNetwork(final IScope scope) throws GamaRuntimeException {
		getCurrentSimulation(scope);
		if(currentSimulation.listNetwork == null){
			currentSimulation.listNetwork = new HashMap<String, MovingOnNetworkSkill.DataSimulation.DataNetwork>();
		}

		DataSimulation.DataNetwork dn = currentSimulation.new DataNetwork();
		currentSimulation.currentDataNetwork = dn;
		dn.lengthAttribute = (String) scope.getArg(IKeywordMoNAdditional.LENGTH_ATTRIBUTE, IType.STRING);;
		dn.speedAttribute = (String) scope.getArg(IKeywordMoNAdditional.SPEED_ATTRIBUTE, IType.STRING);;
		String name = (String) scope.getArg(IKeywordMoNAdditional.NAME, IType.STRING);
		GamaGraph gamaGraph = (GamaGraph) scope.getArg(IKeywordMoNAdditional.GRAPH, IType.GRAPH);
		dn.gamaGraph = gamaGraph;
		Graph graph = new MultiGraph(name, true, false);
		dn.graph = graph;
		Dijkstra dijkstra = null;
		dn.dijkstra = dijkstra;
		FileSinkDGSFiltered fileSink = new FileSinkDGSFiltered();
		dn.fileSink = fileSink;
		//graph.addSink(fileSink);
		String fileName = "";
		try {
			fileName = "./Network_"+name+".dgs";
			File yourFile = new File(fileName);
			yourFile.createNewFile(); // if file already exists will do nothing
		} catch (IOException e) {
			e.printStackTrace();
		}

		try {
			fileSink.begin(fileName);
		} catch (IOException e) {
			e.printStackTrace();
		}

		// Filter useless attributes in edges
		fileSink.addEdgeAttributeFiltered("gama_agent");
		fileSink.addEdgeAttributeFiltered("color");
		fileSink.addEdgeAttributeFiltered("graphstream_edge");
		fileSink.addEdgeAttributeFiltered("gama_time");

		// No need to save graph attributes
		fileSink.setNoFilterGraphAttributeAdded(false);
		fileSink.setNoFilterGraphAttributeChanged(false);
		fileSink.setNoFilterGraphAttributeRemoved(false);

		// and no need either of result which contains Dijsktra reference
		fileSink.addNodeAttributeFiltered("result");
		graph.stepBegins(currentSimulation.currentCycle);
		getGraphstreamGraphFromGamaGraph(scope, currentSimulation, gamaGraph, graph);

		try {
			fileSink.flush();
		} catch (IOException e) {
			e.printStackTrace();
		}
		currentSimulation.listNetwork.put(name, dn);
	}

	@action(
			name = "compute_path_length2",
			args = {
					@arg(name = IKeywordMoNAdditional.SOURCE, type = IType.GEOMETRY , optional = false, doc = @doc("the source of the path.")),
					@arg(name = IKeywordMoNAdditional.TARGET, type = IType.GEOMETRY , optional = false, doc = @doc("the target of the path."))
			},
			doc =
			@doc(value = "Compute the path length", examples = { @example("") })
			)
	public double computeLengthPathAction(final IScope scope) throws GamaRuntimeException {
		final ILocation source = getSource(scope);
		final ILocation target = getTarget(scope);
		return computeShortestPath(scope, source, target);
	}

	@action(
			name = "go_to",
			args = {
					@arg(name = IKeywordMoNAdditional.TARGET, type = IType.GEOMETRY , optional = false, doc = @doc("the location or entity towards which to move.")),
					@arg(name = IKeywordMoNAdditional.MARK, type = IType.FLOAT, optional = true, doc = @doc("The mark (a value) left on the agent's route.")),
			},
			doc =
			@doc(value = "moves the agent towards the target passed in the arguments.", returns = "the path followed by the agent.", examples = { @example("do goto target: (one_of road).location on: road_network;") })
			)
	public double gotoAction(final IScope scope) throws GamaRuntimeException {
		getCurrentSimulation(scope);
		final IAgent agent = getCurrentAgent(scope);
		init(scope, agent);

		// The source is the current location of the current agent
		final ILocation source = agent.getLocation().copy(scope);
		// The target is the location of the thing passing through argument (an agent or a point or a geometry)
		final ILocation target = getTarget(scope);

		if(currentSimulation.currentTarget == null || !currentSimulation.currentTarget.equals(target)){
			// Need to compute the path
			agent.setAttribute("pathLength", computeShortestPath(scope, source, target));
			currentSimulation.currentTarget = target;
			currentSimulation.agentFromOutsideToInside = true;
			currentSimulation.agentInside = false;
			currentSimulation.agentFromInsideToOutside = false;
			currentSimulation.agentOnANode = false;
			currentSimulation.indexSegment = -1;
		}
		// The path has been computed, we need to know how many time the agent has in order to make the move.
		currentSimulation.remainingTime = scope.getClock().getStepInSeconds(); // TODO : check if the value of step should be in second

		if(currentSimulation.currentCycle != scope.getClock().getCycle()){
			currentSimulation.currentCycle = scope.getClock().getCycle();
			for(DataSimulation.DataNetwork dn : currentSimulation.listNetwork.values()){
				// Color the Gama network according to the flow let by agents

				// ===============================================
				// TODO the way to color the graph (which variable is considered and the 0.95%) should be parameters of go_to action available from GAMA for the user
				// ===============================================

				colorGamaGraph(dn.graph, "current_marks");// Available arguments : //cumulative_nb_agents//current_marks//current_nb_agents//cumulative_marks
				dn.graph.stepBegins(currentSimulation.currentCycle);
				for(Edge e : dn.graph.getEachEdge()){
					e.setAttribute("current_volume", 0);
					((IAgent)(e.getAttribute("gama_agent"))).setAttribute("current_volume", e.getNumber("current_volume"));
					if(e.getNumber("current_marks") != 0) {
						e.setAttribute("current_marks", e.getNumber("current_marks")*0.95);
						((IAgent)(e.getAttribute("gama_agent"))).setAttribute("current_marks", e.getNumber("current_marks"));
					}
					if(e.getNumber("current_nb_agents") != 0) {
						e.setAttribute("current_nb_agents", e.getNumber("current_nb_agents")*0.95);
						((IAgent)(e.getAttribute("gama_agent"))).setAttribute("current_nb_agents", e.getNumber("current_nb_agents"));
					}
				}
			}
		}

		IList gl;
		if(currentSimulation.currentGsPathEdge.size()==0 && currentSimulation.agentFromOutsideToInside){
			reachAndLeave(scope, agent, target);
			gl = null;
		}
		else{
			// Move the agent from outside the network to inside (when he is not already on the network).
			movingFromOutsideToInside(scope, agent);

			// Move the agent inside the network (and get the agent edges that the agent has traveled) (when the agent is already on the network and can still move).
			gl = movingInside(scope, agent, target);

			// Move the agent from inside the network to outside (when the target must and can leave the network).
			movingFromInsideToOutside(scope, agent, target);
		}

		agent.setAttribute("agentFromOutsideToInside", currentSimulation.agentFromOutsideToInside);
		agent.setAttribute("agentInside", currentSimulation.agentInside);
		agent.setAttribute("agentFromInsideToOutside", currentSimulation.agentFromInsideToOutside);
		agent.setAttribute("agentOnANode", currentSimulation.agentOnANode);
		agent.setAttribute("indexSegment", currentSimulation.indexSegment);
		agent.setAttribute("currentGsPathEdge", currentSimulation.currentGsPathEdge);
		agent.setAttribute("currentGsPathNode", currentSimulation.currentGsPathNode);
		agent.setAttribute("currentTarget", currentSimulation.currentTarget);

		try {
			currentSimulation.currentDataNetwork.fileSink.flush();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return 0.0;
	}

	private void init(final IScope scope, final IAgent agent) {
		currentSimulation.agentFromOutsideToInside = true;
		if(agent.hasAttribute("agentFromOutsideToInside"))
			currentSimulation.agentFromOutsideToInside = (Boolean) agent.getAttribute("agentFromOutsideToInside");

		currentSimulation.agentInside = false;
		if(agent.hasAttribute("agentInside"))
			currentSimulation.agentInside = (Boolean) agent.getAttribute("agentInside");

		currentSimulation.agentFromInsideToOutside = false;
		if(agent.hasAttribute("agentFromInsideToOutside"))
			currentSimulation.agentFromInsideToOutside = (Boolean) agent.getAttribute("agentFromInsideToOutside");

		currentSimulation.agentOnANode = false;
		if(agent.hasAttribute("agentOnANode"))
			currentSimulation.agentOnANode = (Boolean) agent.getAttribute("agentOnANode");

		currentSimulation.indexSegment = -1;
		if(agent.hasAttribute("indexSegment"))
			currentSimulation.indexSegment = (Integer) agent.getAttribute("indexSegment");

		currentSimulation.currentGsPathEdge = null;
		if(agent.hasAttribute("currentGsPathEdge"))
			currentSimulation.currentGsPathEdge = (List<Edge>) agent.getAttribute("currentGsPathEdge");

		currentSimulation.currentGsPathNode = null;
		if(agent.hasAttribute("currentGsPathNode"))
			currentSimulation.currentGsPathNode = (List<Node>) agent.getAttribute("currentGsPathNode");

		currentSimulation.currentTarget = null;// We use this variable to know if we already have computed the shortest path
		if(agent.hasAttribute("currentTarget"))
			currentSimulation.currentTarget = (ILocation) agent.getAttribute("currentTarget");

		String name = (String)agent.getAttribute("networkType");
		currentSimulation.currentDataNetwork = currentSimulation.listNetwork.get(name);
	}

	private void reachAndLeave(final IScope scope, final IAgent agent, final ILocation target){
		if(currentSimulation.remainingTime > 0){
			Coordinate dest;
			if(!currentSimulation.agentInside){
				dest = new Coordinate(currentSimulation.currentGsPathNode.get(0).getNumber("x"), currentSimulation.currentGsPathNode.get(0).getNumber("y"));
			}
			else{
				dest = new Coordinate(target.getX(), target.getY());
			}

			// The position of the agent
			GamaPoint currentLocation = (GamaPoint) agent.getLocation().copy(scope);
			double xc = currentLocation.getX();
			double yc = currentLocation.getY();

			double dist = Math.hypot(xc - dest.x, yc - dest.y);
			double time = dist / getDefaultSpeed(agent);// We use the default speed because the agent is not yet on the network

			if(currentSimulation.remainingTime >= time){
				currentLocation.setLocation(dest.x, dest.y, currentLocation.getZ());
				agent.setLocation(currentLocation);
				currentSimulation.agentInside = true;
			}
			else{
				double coef = currentSimulation.remainingTime / time;
				double x_inter = xc + (dest.x - xc) * coef;
				double y_inter = yc + (dest.y - yc) * coef;
				currentLocation.setLocation(x_inter, y_inter, currentLocation.getZ());
				agent.setLocation(currentLocation);
			}

			currentSimulation.remainingTime -= time;
		}
	}

	private void movingFromOutsideToInside(final IScope scope, final IAgent agent){
		if(currentSimulation.agentFromOutsideToInside && currentSimulation.remainingTime > 0){
			/*
			 *  First step : find the closest segment to the agent
			 *  Indeed, one edge of the path could be made with more than one segment
			 */
			// The position of the agent
			GamaPoint currentLocation = (GamaPoint) agent.getLocation().copy(scope);
			Point currentPointLocation = (Point) agent.getLocation().copy(scope).getInnerGeometry();
			// The closest road
			IAgent gamaRoad = currentSimulation.currentGsPathEdge.get(0).getAttribute("gama_agent");

			// Find the closest segment among the road's
			double distAgentToNetwork = Double.MAX_VALUE;
			Coordinate coords[] = gamaRoad.getInnerGeometry().getCoordinates();
			Coordinate[] tempCoord = new Coordinate[2];
			int indexBestSegment = 0;
			for ( int i = 0; i < coords.length - 1; i++ ) {
				tempCoord[0] = coords[i];
				tempCoord[1] = coords[i + 1];
				LineString segment = GeometryUtils.GEOMETRY_FACTORY.createLineString(tempCoord);
				double distS = segment.distance(currentPointLocation);
				if ( distS < distAgentToNetwork ) {
					distAgentToNetwork = distS;
					indexBestSegment = i;
				}
			}

			/*
			 * Second step : Find the closest point on this segment
			 */
			// Get coordinates of these different points
			Coordinate dest = getClosestLocation(new Coordinate(currentPointLocation.getX(), currentPointLocation.getY()), coords[indexBestSegment], coords[indexBestSegment+1]);
			double xc = currentPointLocation.getX();
			double yc = currentPointLocation.getY();

			/*
			 * Third step : move the agent on this point
			 */
			double dist = Math.hypot(xc - dest.x, yc - dest.y);
			double time = dist / getDefaultSpeed(agent);// We use the default speed because the agent is not yet on the network

			if(currentSimulation.remainingTime >= time){
				currentLocation.setLocation(dest.x, dest.y, currentLocation.getZ());
				agent.setLocation(currentLocation);
				currentSimulation.agentFromOutsideToInside = false;
				currentSimulation.agentInside = true;
			}
			else{
				double coef = currentSimulation.remainingTime / time;
				double x_inter = xc + (dest.x - xc) * coef;
				double y_inter = yc + (dest.y - yc) * coef;
				currentLocation.setLocation(x_inter, y_inter, currentLocation.getZ());
				agent.setLocation(currentLocation);
			}

			currentSimulation.remainingTime -= time;
		}
	}

	private IList movingInside(final IScope scope, final IAgent agent, final ILocation target){
		if(currentSimulation.agentInside && currentSimulation.remainingTime > 0){
			double mark = (Double) scope.getArg(IKeywordMoNAdditional.MARK, IType.FLOAT);
			ILocation currentLocation = (ILocation) agent.getLocation().copy(scope);
			// It follows the path on the graph, node by node
			IList gl = GamaListFactory.create();

			// Does the agent need to reach the next Node?
			if(!currentSimulation.agentOnANode){
				moveAlongEdge(scope, agent, target, currentSimulation.currentGsPathEdge.get(0));
				// Moreover, if the agent is at the last part of its path, and if he has some remaining time, then, it means that he will leave the network
				if(currentSimulation.currentGsPathEdge.size()== 1 && currentSimulation.remainingTime >= 0){
					// Thus, we pop the current edge of the path and the node (the last ones)
					gl.addValue(scope, currentSimulation.currentGsPathEdge.remove(0).getAttribute("gama_agent"));
					currentSimulation.currentGsPathNode.remove(0);
					currentSimulation.agentFromInsideToOutside = true;
					currentSimulation.agentInside = false;
				}
			}

			if(currentSimulation.currentGsPathNode.isEmpty()){
				currentSimulation.agentFromInsideToOutside = true;
				currentSimulation.agentInside = false;
			}

			while(currentSimulation.remainingTime > 0 && !currentSimulation.currentGsPathEdge.isEmpty()){
				Edge edge = currentSimulation.currentGsPathEdge.get(0);
				double time = edge.getNumber(currentSimulation.currentDataNetwork.lengthAttribute)*100000 / (edge.getNumber(currentSimulation.currentDataNetwork.speedAttribute)*1000/3600);
				// currentGsPath.size()== 1 when the agent is at the end of the path,
				// therefore it must stop before the next node in order to leave the network
				// But he can also stop before the end of an edge if he has not enough remaining time
				if(currentSimulation.currentGsPathEdge.size()== 1 || currentSimulation.remainingTime < time){
					// The moving agent stops between two nodes somewhere on the edge
					// Move the agent to this "somewhere"
					moveAlongEdge(scope, agent, target, edge);
					// Moreover, if the agent is at the last part of its path, and if he has some remaining time, then, it means that he will leave the network
					if(currentSimulation.currentGsPathEdge.size()== 1 && currentSimulation.remainingTime >= 0){
						// Thus, we pop the current edge of the path and the node (the last ones)
						currentSimulation.currentGsPathEdge.remove(0);
						currentSimulation.currentGsPathNode.remove(0);
						currentSimulation.agentFromInsideToOutside = true;
						currentSimulation.agentInside = false;
					}
				}
				else{
					// We move the agent to the next node
					currentSimulation.agentOnANode = true;
					currentSimulation.currentGsPathEdge.get(0).setAttribute("cumulative_marks", currentSimulation.currentGsPathEdge.get(0).getNumber("cumulative_marks") + mark);
					((IAgent)(currentSimulation.currentGsPathEdge.get(0).getAttribute("gama_agent"))).setAttribute("cumulative_marks", currentSimulation.currentGsPathEdge.get(0).getNumber("cumulative_marks"));
					currentSimulation.currentGsPathEdge.get(0).setAttribute("current_marks", currentSimulation.currentGsPathEdge.get(0).getNumber("current_marks") + mark);
					((IAgent)(currentSimulation.currentGsPathEdge.get(0).getAttribute("gama_agent"))).setAttribute("current_marks", currentSimulation.currentGsPathEdge.get(0).getNumber("current_marks"));
					currentSimulation.currentGsPathEdge.get(0).setAttribute("current_volume", currentSimulation.currentGsPathEdge.get(0).getNumber("current_volume") + mark);
					((IAgent)(currentSimulation.currentGsPathEdge.get(0).getAttribute("gama_agent"))).setAttribute("current_volume", currentSimulation.currentGsPathEdge.get(0).getNumber("current_volume"));
					currentSimulation.currentGsPathEdge.get(0).setAttribute("cumulative_nb_agents", currentSimulation.currentGsPathEdge.get(0).getNumber("cumulative_nb_agents") + 1);
					((IAgent)(currentSimulation.currentGsPathEdge.get(0).getAttribute("gama_agent"))).setAttribute("cumulative_nb_agents", currentSimulation.currentGsPathEdge.get(0).getNumber("cumulative_nb_agents"));
					currentSimulation.currentGsPathEdge.get(0).setAttribute("current_nb_agents", currentSimulation.currentGsPathEdge.get(0).getNumber("current_nb_agents") + 1);
					((IAgent)(currentSimulation.currentGsPathEdge.get(0).getAttribute("gama_agent"))).setAttribute("current_nb_agents", currentSimulation.currentGsPathEdge.get(0).getNumber("current_nb_agents"));
					currentSimulation.currentGsPathEdge.remove(0);
					// Set the location of the agent to the next node
					if(currentSimulation.currentGsPathNode.get(0).hasAttribute("gama_agent"))
						currentLocation = ((IAgent)currentSimulation.currentGsPathNode.get(0).getAttribute("gama_agent")).getLocation().toGamaPoint();
					else
						currentLocation = new GamaPoint( currentSimulation.currentGsPathNode.get(0).getNumber("x"), currentSimulation.currentGsPathNode.get(0).getNumber("y"));
					currentSimulation.currentGsPathNode.remove(0);
					//We set the location of the agent in order to make the move
					agent.setLocation(currentLocation);
					currentSimulation.remainingTime -= time;
				}
				// We add the gama agent associated to this edge
				gl.addValue(scope, edge.getAttribute("gama_agent"));
			}
			// We return the list of edges that the agent has traveled.
			return gl;
		}
		// The agent can't move within the network. We return an empty list
		return null;
	}

	private void movingFromInsideToOutside(final IScope scope, final IAgent agent, final ILocation target){
		if(currentSimulation.agentFromInsideToOutside && currentSimulation.remainingTime > 0){
			GamaPoint currentLocation = (GamaPoint) agent.getLocation().copy(scope);
			// Compute the time needed to go to the next side of the segment
			double x_agent = currentLocation.getX();
			double y_agent = currentLocation.getY();

			double x_target = target.getX();
			double y_target = target.getY();

			if(x_agent != x_target || y_agent != y_target){
				double dist = Math.hypot(x_agent - x_target, y_agent - y_target);
				double time = dist / getDefaultSpeed(agent);

				// Move the agent
				if(currentSimulation.remainingTime >= time){
					currentLocation.setLocation(x_target, y_target, currentLocation.getZ());
					currentSimulation.currentTarget = null;
				}
				else{
					double coef = currentSimulation.remainingTime / time;
					double x_inter = x_agent + (x_target - x_agent) * coef;
					double y_inter = y_agent + (y_target - y_agent) * coef;
					currentLocation.setLocation(x_inter, y_inter, currentLocation.getZ());
				}
				agent.setLocation(currentLocation);
				currentSimulation.remainingTime -= time;
			}
		}
	}

	/**
	 * Move the agent along an edge. There are two possibilities to stop the agent on this edges :
	 * - firstly, there is not enough remaining time to reach the end of this edge.
	 * - secondly, the next move will be to leave the network and reach the target (it is also possible that the agent has not enough time to reach the exit point).
	 * @param scope the current scope
	 * @param agent the current agent
	 * @param target the final target to reach
	 * @param e the current edge
	 */
	private void moveAlongEdge(final IScope scope, final IAgent agent, final ILocation target, Edge e){
		double mark = (Double) scope.getArg(IKeywordMoNAdditional.MARK, IType.FLOAT);
		e.setAttribute("cumulative_marks", e.getNumber("cumulative_marks") + mark);
		((IAgent)(e.getAttribute("gama_agent"))).setAttribute("cumulative_marks", e.getNumber("cumulative_marks") + mark);
		e.setAttribute("current_marks", e.getNumber("current_marks"));
		((IAgent)(e.getAttribute("gama_agent"))).setAttribute("current_marks", e.getNumber("current_marks") + mark);
		e.setAttribute("current_volume", e.getNumber("current_volume"));
		((IAgent)(e.getAttribute("gama_agent"))).setAttribute("current_volume", e.getNumber("current_volume") + mark);
		e.setAttribute("cumulative_nb_agents", e.getNumber("cumulative_nb_agents"));
		((IAgent)(e.getAttribute("gama_agent"))).setAttribute("cumulative_nb_agents", e.getNumber("cumulative_nb_agents") + 1);
		e.setAttribute("current_nb_agents", e.getNumber("current_nb_agents") + 1);
		((IAgent)(e.getAttribute("gama_agent"))).setAttribute("current_nb_agents", e.getNumber("current_nb_agents"));

		GamaPoint currentLocation = (GamaPoint) agent.getLocation().copy(scope);
		currentSimulation.agentOnANode = false;
		// Get the geometry of the edge
		IShape shape = ((IAgent)(e.getAttribute("gama_agent"))).getGeometry();
		final Coordinate coords[] = shape.getInnerGeometry().getCoordinates();

		// Determine in which way we must browse the list of segments
		boolean incrementWay = true;
		if(coords[coords.length-1].x != currentSimulation.currentGsPathNode.get(0).getNumber("x") || coords[coords.length-1].y != currentSimulation.currentGsPathNode.get(0).getNumber("y")){
			incrementWay = false;
			if(currentSimulation.indexSegment == -1)
				currentSimulation.indexSegment = coords.length-2;
		}
		else{
			if(currentSimulation.indexSegment == -1)
				currentSimulation.indexSegment = 1;
		}
		// Determine if the agent must stop before the end of the edge because he must leave the network
		// If yes, we look for the segment that the agent must leave
		int indexClosestSegmentToTarget = -1;
		if(e.equals(currentSimulation.currentGsPathEdge.get(currentSimulation.currentGsPathEdge.size()-1))){
			double distTargetToNetwork = Double.MAX_VALUE;
			Coordinate[] tempCoord = new Coordinate[2];
			for ( int j = 0; j < coords.length - 1; j++ ) {
				tempCoord[0] = coords[j];
				tempCoord[1] = coords[j + 1];
				LineString segment = GeometryUtils.GEOMETRY_FACTORY.createLineString(tempCoord);
				double distS = segment.distance(target.getInnerGeometry());
				if ( distS < distTargetToNetwork ) {
					distTargetToNetwork = distS;
					indexClosestSegmentToTarget = j;
				}
			}
		}
		// Browse the segment and move progressively the agent on the edge
		while(currentSimulation.remainingTime >= 0 && currentSimulation.indexSegment >= 0 && currentSimulation.indexSegment < coords.length){
			Coordinate dest;
			if(indexClosestSegmentToTarget != -1 && (currentSimulation.indexSegment == indexClosestSegmentToTarget || (currentSimulation.indexSegment+1) == indexClosestSegmentToTarget) )
				dest = getClosestLocation(new Coordinate(target.getX(), target.getY()), coords[currentSimulation.indexSegment], coords[currentSimulation.indexSegment+1]);
			else 
				dest = coords[currentSimulation.indexSegment];

			// Compute the time needed to go to the next side of the segment
			double x_agent = currentLocation.getX();
			double y_agent = currentLocation.getY();

			double dist = Math.sqrt((x_agent - dest.x)*(x_agent - dest.x) + (y_agent - dest.y)*(y_agent - dest.y));
			double time = dist / (e.getNumber(currentSimulation.currentDataNetwork.speedAttribute)*1000/3600);

			// Move the agent
			if(currentSimulation.remainingTime >= time){
				currentLocation.setLocation(dest.x, dest.y, currentLocation.getZ());
				// Increment or decrement i according to the way that we browse the list of segments
				if(incrementWay)
					currentSimulation.indexSegment++;
				else
					currentSimulation.indexSegment--;
			}
			else{
				double coef = currentSimulation.remainingTime/time;
				double x_inter = x_agent + (dest.x-x_agent)*coef ;
				double y_inter = y_agent + (dest.y-y_agent)*coef ;
				currentLocation.setLocation(x_inter, y_inter, currentLocation.getZ());
			}
			agent.setLocation(currentLocation);
			// Update the remaining time
			currentSimulation.remainingTime -= time;
		}

		if(agent.getLocation().getX() == currentSimulation.currentGsPathNode.get(0).getNumber("x") && agent.getLocation().getY() == currentSimulation.currentGsPathNode.get(0).getNumber("y")){
			currentSimulation.agentOnANode = true;
			currentSimulation.currentGsPathEdge.remove(0);
			currentSimulation.currentGsPathNode.remove(0);
			currentSimulation.indexSegment = -1;
		}
	}

	private Coordinate getClosestLocation(Coordinate coordOutNetwork, Coordinate a, Coordinate b){
		// Get coordinates of these different points
		double xa = a.x;
		double ya = a.y;
		double xb = b.x;
		double yb = b.y;
		double xc = coordOutNetwork.x;
		double yc = coordOutNetwork.y;

		// Compute coordinates of vectors
		// CA Vector
		double ACy = yc - ya;
		double ACx = xc - xa;
		// AB vector
		double ABy = yb - ya;
		double ABx = xb - xa;
		// CB vector
		double BCy = yc - yb;
		double BCx = xc - xb;
		// BA vector
		double BAy = ya - yb;
		double BAx = xa - xb;

		// Compute the angles
		// The angle between ->AC and ->AB
		double CAB = Math.abs( Math.toDegrees(Math.atan2(ACy, ACx)-Math.atan2(ABy, ABx)) );
		// The angle between ->BC and ->BA
		double CBA = Math.abs( Math.toDegrees(Math.atan2(BCy, BCx)-Math.atan2(BAy, BAx)) );

		// Let A and B the nodes of this segment and C be the currentLocation
		// If one of the angles CAB or CBA  is obtuse ( ie.  90 < CAB < 180 or 90 < CBA < 180)
		// 	then the next location is on the segment between C and A (or C and B)
		double x_dest;
		double y_dest;
		if(CAB >= 90 ){
			// Between C and A
			x_dest = xa;
			y_dest = ya;
		}
		else if(CBA >= 90){
			// Between C and B
			x_dest = xb;
			y_dest = yb;
		}
		else {
			// Let H be the orthographic projection of C on AB (thus we have : (CH) _|_ (AB) )
			// The next location is on the segment between C and H
			// Compute unit vector
			double xv = (xb-xa);
			double yv = (yb- ya);
			// Compute distance
			double AH = ( (xc-xa)*xv + (yc-ya)*yv ) / ( Math.sqrt(xv*xv +yv*yv) );
			x_dest = xa + ( AH / (Math.sqrt(xv*xv +yv*yv)) ) * xv;
			y_dest = ya + ( AH / (Math.sqrt(xv*xv +yv*yv)) ) * yv;
		}

		return new Coordinate(x_dest, y_dest);
	}

	private double computeShortestPath(final IScope scope, ILocation source, ILocation target){
		if(currentSimulation.currentDataNetwork.dijkstra == null){
			currentSimulation.currentDataNetwork.dijkstra = new Dijkstra(Dijkstra.Element.EDGE, "result", "gama_time");
			currentSimulation.currentDataNetwork.dijkstra.init(currentSimulation.currentDataNetwork.graph);
		}

		/*
		 *  Find the graphstream source and target node
		 */
		GraphTopology gt = (GraphTopology)(Cast.asTopology(scope, currentSimulation.currentDataNetwork.gamaGraph));
		ITopology topo = scope.getSimulation().getTopology();
		// Find the source node
		IAgent gamaSourceEdge = topo.getAgentClosestTo(scope, source, In.edgesOf(gt.getPlaces()));//gt.getAgentClosestTo(scope, source, In.edgesOf(gt.getPlaces()));
		Edge gsSourceEdge = (Edge)gamaSourceEdge.getAttribute("mon_graphstream_edge");
		Node sourceNode = gsSourceEdge.getNode0();
		// Find the target node
		IAgent gamaTargetEdge = topo.getAgentClosestTo(scope, target, In.edgesOf(gt.getPlaces()));//gt.getAgentClosestTo(scope, target, In.edgesOf(gt.getPlaces()));
		Edge gsTargetEdge = (Edge)gamaTargetEdge.getAttribute("mon_graphstream_edge");
		Node targetNode = gsTargetEdge.getNode0();

		/*
		 *  Compute and get the path
		 */
		currentSimulation.currentDataNetwork.dijkstra.setSource(sourceNode);
		currentSimulation.currentDataNetwork.dijkstra.compute();
		Path p = currentSimulation.currentDataNetwork.dijkstra.getPath(targetNode);

		double length = p.getPathWeight("gama_time");
		final IAgent agent = getCurrentAgent(scope);
		agent.setAttribute(IKeywordMoNAdditional.PATH_LENGTH, length);

		/*
		 * Add closest edge(s)
		 */
		// Add closest source edge to the path if it is missing
		if(!p.contains(gsSourceEdge)){
			sourceNode = gsSourceEdge.getNode1();
			currentSimulation.currentDataNetwork.dijkstra.setSource(sourceNode);
			currentSimulation.currentDataNetwork.dijkstra.compute();
			p = currentSimulation.currentDataNetwork.dijkstra.getPath(targetNode);
		}

		// Add closest target edge to the path if it is missing
		if(!p.contains(gsTargetEdge)){
			targetNode = gsTargetEdge.getNode1();
			p = currentSimulation.currentDataNetwork.dijkstra.getPath(targetNode);
		}

		currentSimulation.currentGsPathEdge = p.getEdgePath();
		currentSimulation.currentGsPathNode = p.getNodePath();
		if(currentSimulation.currentGsPathEdge.size() != 0)
			currentSimulation.currentGsPathNode.remove(0);// The first node is useless

		return p.getPathWeight("gama_time");
	}

	/**
	 * Takes a gama graph as an input, returns a graphstream graph as
	 * close as possible. Preserves double links (multi graph).
	 * Copy of the method of GraphUtilsGraphStream but we save the gama agent in each edges/nodes and the graphstream edge in each gama edge agent
	 * @param gamaGraph
	 * @return The Graphstream graph
	 */
	private static void getGraphstreamGraphFromGamaGraph(IScope scope, DataSimulation currentSimulation, final IGraph gamaGraph, Graph g) {
		Map<Object, Node> gamaNode2graphStreamNode = new HashMap<Object, Node>(gamaGraph._internalNodesSet().size());
		// add nodes
		for ( Object v : gamaGraph._internalVertexMap().keySet() ) {
			_Vertex vertex = (_Vertex) gamaGraph._internalVertexMap().get(v);
			Node n = g.addNode(v.toString());
			gamaNode2graphStreamNode.put(v, n);
			if ( v instanceof IAgent ) {
				IAgent a = (IAgent) v;
				n.addAttribute("gama_agent", a);
				for ( Object key : a.getOrCreateAttributes().keySet() ) {
					Object value = GraphUtilsGraphStream.preprocessGamaValue(a.getOrCreateAttributes().get(key));
					if(value != null)
						n.addAttribute(key.toString(), value.toString());
				}
			}

			if ( v instanceof IShape ) {
				IShape sh = (IShape) v;
				n.setAttribute("x", sh.getLocation().getX());
				n.setAttribute("y", sh.getLocation().getY());
				n.setAttribute("z", sh.getLocation().getZ());
			}
		}
		// add edges
		for ( Object edgeObj : gamaGraph._internalEdgeMap().keySet() ) {
			_Edge edge = (_Edge) gamaGraph._internalEdgeMap().get(edgeObj);
			try {
				Edge e = // We call the function where we give the nodes object directly, is it more efficient than give the string id? Because, if no, we don't need the "gamaNode2graphStreamNode" map...
						g.addEdge(edgeObj.toString(), gamaNode2graphStreamNode.get(edge.getSource()), gamaNode2graphStreamNode.get(edge.getTarget()),
								gamaGraph.isDirected() );// till now, directionality of an edge depends on the whole gama graph
				if ( edgeObj instanceof IAgent ) {
					IAgent a = (IAgent) edgeObj;
					// e know a
					e.addAttribute("gama_agent", a);
					for ( Object key : a.getOrCreateAttributes().keySet() ) {
						Object value = GraphUtilsGraphStream.preprocessGamaValue(a.getOrCreateAttributes().get(key));
						if(value != null)
							e.addAttribute(key.toString(), value.toString());
					}
					e.addAttribute("gama_time", e.getNumber(currentSimulation.currentDataNetwork.lengthAttribute) * e.getNumber(currentSimulation.currentDataNetwork.speedAttribute));
					e.setAttribute("current_marks", 0.0);
					e.setAttribute("cumulative_marks", 0.0);
					e.setAttribute("cumulative_nb_agents", 0.0);
					e.setAttribute("current_nb_agents", 0.0);
					// a know e
					a.setAttribute("mon_graphstream_edge", e);
				}
			} catch (EdgeRejectedException e) {
				GAMA.reportError(scope, GamaRuntimeException
						.warning("an edge was rejected during the transformation, probably because it was a double one", scope),
						true);
			} catch (IdAlreadyInUseException e) {
				GAMA.reportError(scope, GamaRuntimeException
						.warning("an edge was rejected during the transformation, probably because it was a double one", scope),
						true);
			}

		}
		// some basic tests for integrity
		if ( gamaGraph.getVertices().size() != g.getNodeCount() ) {
			GAMA.reportError(scope,
					GamaRuntimeException.warning("The exportation ran without error, but an integrity test failed: " +
							"the number of vertices is not correct(" + g.getNodeCount() + " instead of " +
							gamaGraph.getVertices().size() + ")", scope), true);
		}
		if ( gamaGraph.getEdges().size() != g.getEdgeCount() ) {
			GAMA.reportError(scope,
					GamaRuntimeException.warning("The exportation ran without error, but an integrity test failed: " +
							"the number of edges is not correct(" + g.getEdgeCount() + " instead of " +
							gamaGraph.getEdges().size() + ")", scope), true);
		}
	}

	/**
	 * Color the Gama graph according to the given attributes' value on the edges
	 * @param attr The attribute containing the value allowing the coloring
	 */
	private void colorGamaGraph(Graph g, String attr){
		SortedList listEdge = new SortedList();
		for(Edge e: g.getEachEdge()) {
			listEdge.add(e, attr);
		}
		listEdge.sort();
		for(int i = 0; i < listEdge.size(); i++){
			for(Edge e : listEdge.get(i).getEdges() ){
				if(listEdge.size() != 1){
					// Color edges
					// 	#fed976 kind of yellow
					// 	#feb24c
					// 	#fd8d3c
					// 	#fc4e2a
					// 	#e31a1c
					// 	#bd0026
					// 	#800026 strong red
					if(i * 10 < (listEdge.size()/7.) * 10) {
						((IAgent)(e.getAttribute("gama_agent"))).setAttribute("colorRVBValue", Integer.parseInt("800026",16));
					}
					else if(i * 10 < (listEdge.size()/7.) * 20) {
						((IAgent)(e.getAttribute("gama_agent"))).setAttribute("colorRVBValue", Integer.parseInt("bd0026",16));
					}
					else if(i * 10 < (listEdge.size()/7.) * 30) {
						((IAgent)(e.getAttribute("gama_agent"))).setAttribute("colorRVBValue", Integer.parseInt("e31a1c",16));
					}
					else if(i * 10 < (listEdge.size()/7.) * 40) {
						((IAgent)(e.getAttribute("gama_agent"))).setAttribute("colorRVBValue", Integer.parseInt("fc4e2a",16));
					}
					else if(i * 10 < (listEdge.size()/7.) * 50) {
						((IAgent)(e.getAttribute("gama_agent"))).setAttribute("colorRVBValue", Integer.parseInt("fd8d3c",16));
					}
					else if(i * 10 < (listEdge.size()/7.) * 60) {
						((IAgent)(e.getAttribute("gama_agent"))).setAttribute("colorRVBValue", Integer.parseInt("feb24c",16));
					}
					else  {
						((IAgent)(e.getAttribute("gama_agent"))).setAttribute("colorRVBValue", Integer.parseInt("fed976",16));
					}
					// And give different size to edges
					if(i * 10 < (listEdge.size()/4.) * 10) {
						((IAgent)(e.getAttribute("gama_agent"))).setAttribute("sizeShape", 3);
					}
					else if(i * 10 < (listEdge.size()/4.) * 20) {
						((IAgent)(e.getAttribute("gama_agent"))).setAttribute("sizeShape", 2);
					}
					else if(i * 10 < (listEdge.size()/4.) * 30) {
						((IAgent)(e.getAttribute("gama_agent"))).setAttribute("sizeShape", 1);
					}
					else  {
						((IAgent)(e.getAttribute("gama_agent"))).setAttribute("sizeShape", 0);
					}
				}
				else {
					((IAgent)(e.getAttribute("gama_agent"))).setAttribute("colorRVBValue", Integer.parseInt("fed976",16));
					((IAgent)(e.getAttribute("gama_agent"))).setAttribute("sizeShape", 1);
				}
			}
		}
	}

	private class SortedList {
		ArrayList<Edges> list;

		public SortedList(){
			list = new ArrayList<MovingOnNetworkSkill.Edges>();
		}

		public boolean add(Edge e, String attr){
			double eVal = 0.0;
			if(e.hasAttribute(attr) && e.getNumber(attr) > 0)
				eVal = e.getNumber(attr);
			for(Edges es : list){
				final double EPSILON = 1;
				if( ((es.value - eVal) * (es.value - eVal) < EPSILON * EPSILON) ){
				//if(es.value - 1 <= eVal && eVal <= es.value + 1){
				//if(eVal == es.value){
					return es.add(e, attr);
				}
			}
			Edges es = new Edges(e, attr);
			return list.add(es);
		}

		public int size(){
			return list.size();
		}

		public Edges get(int i){
			return list.get(i);
		}

		public void sort(){
			Collections.sort(list);
		}
	}

	private class Edges implements Comparable<Edges>{
		double value;
		double sumValues;
		ArrayList<Edge> edges;

		public Edges(Edge e, String attr){
			if(e.hasAttribute(attr) && e.getNumber(attr) > 0)
				this.value = e.getNumber(attr);
			else
				this.value = 0.0;
			this.sumValues = this.value;
			this.edges = new ArrayList<Edge>();
			this.edges.add(e);
		}

		public boolean add(Edge e, String attr){
			if(e.hasAttribute(attr) && e.getNumber(attr) > 0)
				sumValues += e.getNumber(attr);
			value = sumValues / edges.size();
			return edges.add(e);
		}

		public ArrayList<Edge> getEdges(){
			return edges;
		}

		public int compareTo(Edges e) {
			if(e.value == value)
				return 0;
			if(value < e.value)
				return 1;
			else
				return -1;
		}
	}
}