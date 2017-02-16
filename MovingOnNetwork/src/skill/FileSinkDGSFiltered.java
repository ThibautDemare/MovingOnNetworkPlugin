package skill;

import java.awt.Color;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.reflect.Array;
import java.util.HashMap;
import java.util.Locale;

import org.graphstream.graph.CompoundAttribute;

/**
 * <p>
 * File output for the DGS (Dynamic Graph Stream) file format. It includes also the possibility to filter dynamic events such as :
 * <ul>
 * 		<li>
 *   		Addition or deletion of nodes or edges.
 * 		</li>
 * 		<li>
 *   		Addition/deletion/modification of attributes of nodes, edges or the graph itself.
 * 		</li>
 * 		<li>
 *   		A step event.
 * 		</li>
 * </ul>
 * </p>
 *
 * <p>
 * For instance :
 * </p>
 *
 * <pre>
 * Graph graph = new SingleGraph("My_Graph");
 * FileSinkDGSFiltered fileSink = new FileSinkDGSFiltered();
 * graph.addSink(fileSink);
 *
 * // No need to save the attribute "attr1" of any edges
 * fileSink.addEdgeAttributeFiltered("attr1");
 * // No need to save graph attributes
 * fileSink.setNoFilterGraphAttributeAdded(false);
 *
 * // Start to listen event
 * fileSink.begin("./my_graph.dgs");
 *
 * // Make some modifications on the graph and generate events
 * graph.stepBegins(0); // this event will be saved
 * graph.addAttribute("attr2", 2); // this event will not be saved
 * Node a = graph.addNode("A"); // this event will be saved
 * a.addAttribute("attr3", 3); // this event will be saved
 *
 * // and now, no more need to save modification on nodes attributes
 * fileSink.setNoFilterNodeAttributeChanged(false);
 *
 * Node b = graph.addNode("B"); // this event will be saved
 * b.addAttribute("attr4", 4); // this event will not be saved
 *
 * Edge ab = graph.addEdge("AB", a, b); // this event will be saved
 * ab.addAttribute("attr1", 1); // this event will not be saved
 * ab.addAttribute("attr5", 5); // this event will be saved
 *
 * fileSink.end();
 * </pre>
 */
public class FileSinkDGSFiltered extends FileSinkBaseFiltered {

	// Attribute
	
	/**
	 * A shortcut to the output.
	 */
	protected PrintWriter out;
	protected String graphName = "";
	
	//Command
	
	@Override
	protected void outputHeader() throws IOException {
		out = (PrintWriter) output;
		out.printf("DGS004%n");

		if (graphName.length() <= 0)
			out.printf("null 0 0%n");
		else
			out.printf("\"%s\" 0 0%n", formatStringForQuoting(graphName));
	}
	
	@Override
	protected void outputEndOfFile() throws IOException {
		//NOP
	}
	
	public void edgeAttributeAdded(String graphId, long timeId, String edgeId,
			String attribute, Object value) {
		if(noFilterEdgeAttributeAdded && !edgeAttributesFiltered.contains(attribute))
			edgeAttributeChanged(graphId, timeId, edgeId, attribute, null, value);
	}
	
	public void edgeAttributeChanged(String graphId, long timeId,
			String edgeId, String attribute, Object oldValue, Object newValue) {
		if(noFilterEdgeAttributeChanged && !edgeAttributesFiltered.contains(attribute))
			out.printf("ce \"%s\" %s%n", formatStringForQuoting(edgeId),
				attributeString(attribute, newValue, false));
	}
	
	public void edgeAttributeRemoved(String graphId, long timeId,
			String edgeId, String attribute) {
		if(noFilterEdgeAttributeRemoved && !edgeAttributesFiltered.contains(attribute))
			out.printf("ce \"%s\" %s%n", formatStringForQuoting(edgeId),
				attributeString(attribute, null, true));
	}
	
	public void graphAttributeAdded(String graphId, long timeId,
			String attribute, Object value) {
		if(noFilterGraphAttributeAdded && !graphAttributesFiltered.contains(attribute))
			graphAttributeChanged(graphId, timeId, attribute, null, value);
	}
	
	public void graphAttributeChanged(String graphId, long timeId,
			String attribute, Object oldValue, Object newValue) {
		if(noFilterGraphAttributeChanged && !graphAttributesFiltered.contains(attribute))
			out.printf("cg %s%n", attributeString(attribute, newValue, false));
	}
	
	public void graphAttributeRemoved(String graphId, long timeId,
			String attribute) {
		if(noFilterGraphAttributeRemoved && !graphAttributesFiltered.contains(attribute))
			out.printf("cg %s%n", attributeString(attribute, null, true));
	}
	
	public void nodeAttributeAdded(String graphId, long timeId, String nodeId,
			String attribute, Object value) {
		if(noFilterNodeAttributeAdded && !nodeAttributesFiltered.contains(attribute))
			nodeAttributeChanged(graphId, timeId, nodeId, attribute, null, value);
	}
	
	public void nodeAttributeChanged(String graphId, long timeId,
			String nodeId, String attribute, Object oldValue, Object newValue) {
		if(noFilterNodeAttributeChanged && !nodeAttributesFiltered.contains(attribute))
			out.printf("cn \"%s\" %s%n", formatStringForQuoting(nodeId),
				attributeString(attribute, newValue, false));
	}
	
	public void nodeAttributeRemoved(String graphId, long timeId,
			String nodeId, String attribute) {
		if(noFilterNodeAttributeRemoved && !nodeAttributesFiltered.contains(attribute))
			out.printf("cn \"%s\" %s%n", formatStringForQuoting(nodeId),
				attributeString(attribute, null, true));
	}
	
	public void edgeAdded(String graphId, long timeId, String edgeId,
			String fromNodeId, String toNodeId, boolean directed) {
		if(noFilterEdgeAdded){
			edgeId = formatStringForQuoting(edgeId);
			fromNodeId = formatStringForQuoting(fromNodeId);
			toNodeId = formatStringForQuoting(toNodeId);
			out.printf("ae \"%s\" \"%s\" %s \"%s\"%n", edgeId, fromNodeId,
					directed ? ">" : "", toNodeId);
		}
	}
	
	public void edgeRemoved(String graphId, long timeId, String edgeId) {
		if(noFilterEdgeRemoved)
			out.printf("de \"%s\"%n", formatStringForQuoting(edgeId));
	}
	
	public void graphCleared(String graphId, long timeId) {
		if(noFilterGraphCleared)
			out.printf("cl%n");
	}
	
	public void nodeAdded(String graphId, long timeId, String nodeId) {
		if(noFilterNodeAdded)
			out.printf("an \"%s\"%n", formatStringForQuoting(nodeId));
	}
	
	public void nodeRemoved(String graphId, long timeId, String nodeId) {
		if(noFilterNodeRemoved)
			out.printf("dn \"%s\"%n", formatStringForQuoting(nodeId));
	}
	
	public void stepBegins(String graphId, long timeId, double step) {
		if(noFilterStepBegins)
			out.printf(Locale.US, "st %f%n", step);
	}
	
	//Utility
	
	protected String formatStringForQuoting(String str) {
		return str.replaceAll("(^|[^\\\\])\"", "$1\\\\\"");
	}
	
	protected String attributeString(String key, Object value, boolean remove) {
		if (key == null || key.length() == 0)
			return null;
		if (remove) {
			return String.format(" -\"%s\"", key);
		} else {
			if (value != null && value.getClass().isArray())
				return String.format(" \"%s\":%s", key, arrayString(value));
			else
				return String.format(" \"%s\":%s", key, valueString(value));
		}
	}
	
	protected String arrayString(Object value) {
		if (value != null && value.getClass().isArray()) {
			StringBuilder sb = new StringBuilder();
			sb.append("{");
			if (Array.getLength(value) == 0)
				sb.append("\"\"");
			else
				sb.append(arrayString(Array.get(value, 0)));
			for (int i = 1; i < Array.getLength(value); ++i)
				sb.append(String
						.format(",%s", arrayString(Array.get(value, i))));
			sb.append("}");
			return sb.toString();
		} else {
			return valueString(value);
		}
	}
	
	protected String valueString(Object value) {
		if (value == null)
			return "\"\"";
		if (value instanceof CharSequence) {
			if (value instanceof String)
				return String.format("\"%s\"",
						formatStringForQuoting((String) value));
			else
				return String.format("\"%s\"", (CharSequence) value);
		} else if (value instanceof Number) {
			Number nval = (Number) value;
			if (value instanceof Integer || value instanceof Short
					|| value instanceof Byte || value instanceof Long)
				return String.format(Locale.US, "%d", nval.longValue());
			else
				return String.format(Locale.US, "%f", nval.doubleValue());
		} else if (value instanceof Boolean) {
			return String.format(Locale.US, "%b", ((Boolean) value));
		} else if (value instanceof Character) {
			return String.format("\"%c\"", ((Character) value).charValue());
		} else if (value instanceof Object[]) {
			Object array[] = (Object[]) value;
			int n = array.length;
			StringBuffer sb = new StringBuffer();
			if (array.length > 0)
				sb.append(valueString(array[0]));
			for (int i = 1; i < n; i++) {
				sb.append(",");
				sb.append(valueString(array[i]));
			}
			return sb.toString();
		} else if (value instanceof HashMap<?, ?>
		|| value instanceof CompoundAttribute) {
			HashMap<?, ?> hash;
			if (value instanceof CompoundAttribute)
				hash = ((CompoundAttribute) value).toHashMap();
			else
				hash = (HashMap<?, ?>) value;
			return hashToString(hash);
		} else if (value instanceof Color) {
			Color c = (Color) value;
			return String.format("#%02X%02X%02X%02X", c.getRed(), c.getGreen(),
					c.getBlue(), c.getAlpha());
		} else {
			return String.format("\"%s\"", value.toString());
		}
	}
	
	protected String hashToString(HashMap<?, ?> hash) {
		StringBuilder sb = new StringBuilder();
		sb.append("[ ");
		for (Object key : hash.keySet()) {
			sb.append(attributeString(key.toString(), hash.get(key), false));
			sb.append(",");
		}
		sb.append(']');
		return sb.toString();
	}
}
