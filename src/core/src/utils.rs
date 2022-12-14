use indicatif::{ProgressBar, ProgressStyle};

pub fn styled_progress_bar(total: u64, message: &str) -> ProgressBar {
    let progress_bar = ProgressBar::new(total);
    let format_string = format!(
        "[+]\t{}: {{bar:40.cyan/blue}} {{pos:>7}}/{{len:12}}",
        message
    );

    let pb_style_result = ProgressStyle::with_template(format_string.as_str());
    let pb_style = match pb_style_result {
        Ok(style) => style,
        Err(error) => panic!("Problem with the progress bar: {:?}", error),
    };
    progress_bar.set_style(pb_style);
    progress_bar
}
